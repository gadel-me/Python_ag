from __future__ import print_function, division
import os
import copy
import shutil as sl
import re
import argparse
import math
import pdb
import time
import numpy as np
from natsort import natsorted
#import itertools as it
import scipy.stats
from mpi4py import MPI
from lammps import lammps, PyLammps
import Transformations as cgt
import md_elements as mde
import md_box as mdb
import ag_unify_md as agum
import ag_geometry as agm
import ag_lammps as aglmp
import ag_unify_log as agul
import ag_vectalg as agv

#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator

#==============================================================================#
# Helper functions and variables
#==============================================================================#
_sqrt_3 = math.sqrt(3)


class LmpShortcuts(object):
    """
    """
    def __init__(self, tstart=None, tstop=None, pstart=None, pstop=None, logsteps=None, runsteps=None, pc_file=None, settings_file=None, input_lmpdat=None, input_lmprst=None, inter_lmprst=None, output_lmprst=None, output_lmplog=None, output_dcd=None, output_lmpdat=None, output_name=None, gpu=False):
        """
        """
        self.tstart = tstart
        self.tstop = tstop
        self.pstart = pstart
        self.pstop = pstop
        self.logsteps = logsteps
        self.runsteps = runsteps
        self.pc_file = pc_file
        self.settings_file = settings_file
        self.input_lmpdat = input_lmpdat
        self.input_lmprst = input_lmprst
        self.inter_lmprst = inter_lmprst
        self.output_lmplog = output_lmplog
        self.output_dcd = output_dcd
        self.output_lmpdat = output_lmpdat
        self.output_lmprst = output_lmprst
        self.output_name = output_name
        self.gpu = gpu

    def read_system(self, lmp):
        """
        Read a restart file from another stage, previous run or a data file instead.
        """
        # reading order: output restart -> input restart -> data
        if self.output_lmprst is not None and os.path.isfile(self.output_lmprst):
            lmp.command("read_restart {}".format(self.output_lmprst))
        elif self.input_lmprst is not None and os.path.isfile(self.input_lmprst):
            lmp.command("read_restart {}".format(self.input_lmprst))
        else:
            lmp.command("read_data {}".format(self.input_lmpdat))
            lmp.command("velocity all create {} 4928459 mom yes rot yes dist gaussian".format(self.tstart))

    def unfix_undump(self, pylmp, lmp):
        """
        Remove all fixes and dumps.
        """
        lmp_fixes = []
        lmp_dumps = []

        if rank == 0:
            for fix in pylmp.fixes:
                lmp_fixes.append(fix["name"])
            for dump in pylmp.dumps:
                lmp_dumps.append(dump["name"])

        lmp_fixes = comm.bcast(lmp_fixes, 0)
        lmp_dumps = comm.bcast(lmp_dumps, 0)

        for fix in lmp_fixes:
            if "gpu" in fix:
                continue
            lmp.command("unfix {}".format(fix))

        for dump in lmp_dumps:
            lmp.command("undump {}".format(dump))

    def thermo(self, lmp):
        """
        Log thermodynamic data.
        """
        thermargs = ["step", "temp", "press", "vol", "density", "cella", "cellb", "cellc", "cellalpha", "cellbeta", "cellgamma", "etotal", "pe", "evdwl", "ecoul", "ebond", "eangle", "edihed", "eimp", "enthalpy"]
        lmp.command("thermo_style custom " + " ".join(thermargs))
        lmp.command("thermo_modify lost warn flush yes")
        #lmp.command("thermo_modify line multi format float %g")
        lmp.command("thermo {}".format(self.logsteps))

    def dump(self, lmp, unwrap=True):
        """
        """
        # trajectory
        lmp.command("dump trajectory all dcd {} {}".format(self.logsteps, self.output_dcd))
        lmp.command("restart {} {} {}".format(self.logsteps*50, self.inter_lmprst, self.inter_lmprst))

        if unwrap is True:
            lmp.command("dump_modify trajectory unwrap yes")

    def berendsen(self, lmp):
        """
        """
        lmp.command("fix integrator all nve")
        lmp.command("fix thermostat all temp/berendsen {} {} 0.5".format(self.tstart, self.tstop))
        lmp.command("fix barostat all press/berendsen iso {} {} 50".format(self.pstart, self.pstop))

    def use_gpu(self, lmp, neigh=True):
        """
        """
        if neigh:
            lmp.command("package gpu 1 neigh yes")
        else:
            lmp.command("package gpu 1 neigh no")

        lmp.command("suffix gpu")

    def minimize(self, lmp, style="cg", box_relax=False):
        """
        """
        if box_relax:
            lmp.command("fix box_relax all box/relax aniso {}".format(self.pstart))

        lmp.command("min_style {}".format(style))
        lmp.command("min_modify dmax 2.0")
        lmp.command("minimize 1.0e-9 1.0e-12 100000 1000000")


def _molecules_radii(lmp_sys):
    """
    """
    # get center of geometry and sphere of each solvate molecule which is placed
    # inside the solvent box
    radii_sphere = []
    cogs_solvate = []

    for molecule in lmp_sys.molecules:
        sphere, cog = lmp_sys.get_mol_radius(-1, *molecule)
        radii_sphere.append(sphere)
        cogs_solvate.append(cog)

    return (radii_sphere, cogs_solvate)


def mkdir(folder):
    """
    Create folder or skip creation if it already exists.

    Create a folder, if it exists, do nothing.
    """
    try:
        os.mkdir(folder, 0755)
    except OSError:
        print("***Info: Folder {} already exists!".format(folder))


def _read_system(lmpdat, dcd=None, frame_idx=-1):
    """
    Read a lammps data file and a dcd file on top (optional).

    This function simplifies the process of reading a lammps data file and loading
    a dcd file on top of it for further coordinates. It returns a Universe object
    which has many methods to manipulate it with.

    Parameters
    ----------
    lmpdat : str
        Name of the lammps data file

    dcd : str, optional
        Name of the dcd file with the same amounts of atoms as the lammps data
        file

    frame_idx : int, optional
        Index of the frame to use from the dcd file.

    Returns
    -------
    md_sys : LmpStuff object
        An object of LmpStuff which can be further processed.

    """
    # read output from quenching
    md_sys = aglmp.LmpStuff()
    md_sys.read_lmpdat(lmpdat)

    if dcd is not None:
        md_sys.import_dcd(dcd)
        # since we are only interested in one frame, delete all others
        md_sys.ts_coords = []
        md_sys.ts_boxes = []

        # enable reading the last frame with negative indexing
        if frame_idx == -1:
            md_sys.read_frames(frame=frame_idx - 1, to_frame=frame_idx)
        else:
            md_sys.read_frames(frame=frame_idx, to_frame=frame_idx + 1)

        md_sys.close_dcd()

    return md_sys


def merge_sys(sys_a, sys_b, frame_id_a=-1, frame_id_b=-1, pair_coeffs=False):
    """
    """
    sys_both = copy.deepcopy(sys_a)
    sys_both.extend_universe(sys_b, u1_frame_id=frame_id_a, u2_frame_id=frame_id_b, mode="merge")

    if pair_coeffs is True:
        sys_both.mix_pair_types(mode="ii", mix_style="arithmetic")

    sys_both.fetch_molecules_by_bonds()
    sys_both.mols_to_grps()
    return sys_both


def check_aggregate(mdsys, frame_id=-1, atm_atm_dist=4, unwrap=False, debug=False):
    """
    Check if several molecules form an aggregate.

    #TODO Buggy, does not recognize aggregates properly when they form a chain?

    At the moment this works only for orthogonal cells (reason: sub-cell length)
    Check if the aggregate did not get dissolved in the process. This is achieved
    by calculating the center of geometry for each molecule and subsequent com-
    parison of the distance between each center of geometry. The smallest distance
    must not be larger than four times the radius of the biggest molecule in the
    system.

    Input:
        > mdsys             ag_unify_md.Unification; system output from 'SYS-PREPARATION'-step
        > frame_id          int; frame to process
        > atm_atm_dist      float; radius around an atom in which another atom
                            should be positioned
        > debug             boolean; True if further output should be given,
                            default=False

    Output:
        > aggregate_ok      boolean; True, if aggregate is still intact, False if
                            a part of it drifted away
    """
    aggregate_ok = False

    # lammps may internally have the coordinates wrapped -> unwrap cell when
    # necessary (e.g. coordinates were taken from restart)
    if unwrap is True:
        print("***Check Aggregate Info: Unwrapping cell")
        mdsys.unwrap_cell(frame_id)

    # maximal possible distance between two atoms
    cube_side = atm_atm_dist / _sqrt_3  # -> atm_atm_dist/(2*math.sqrt(3)) * 2

    if debug is True:
        print("***Check Aggregate Info: Cube side {}".format(cube_side))

    mdsys.create_linked_cells(frame_id, rcut_a=cube_side, rcut_b=cube_side,
                              rcut_c=cube_side)

    # include same molecule or it gets missing in our aggregates set
    close_atoms, aggregates = mdsys.chk_atm_dist(frame_id,
                                                 min_dist=atm_atm_dist,
                                                 exclude_same_molecule=False,
                                                 get_aggregates=True,
                                                 debug=False)

    if debug is True:
        print("***Check Aggregate Info: Number of cubes in direction a {}".format(mdsys.ts_lnk_cls[frame_id].ra))
        print("***Check Aggregate Info: Number of cubes in direction b {}".format(mdsys.ts_lnk_cls[frame_id].rb))
        print("***Check Aggregate Info: Number of cubes in direction c {}".format(mdsys.ts_lnk_cls[frame_id].rc))
        print("***Group IDs of all aggregates: {}".format(aggregates))
        print("***IDs of close atoms: {}".format(" ".join([str(i) for i in close_atoms])))

    # aggregate is only o.k. if all molecules are part of it
    if len(aggregates) == 1 and len(aggregates[0]) == len(mdsys.molecules):
        aggregate_ok = True

        if debug is True:
            print("***Check Aggregate Info: Aggregate looks fine!")

    return aggregate_ok


def cut_box(lmpdat, box, lmpdat_out, dcd=None):
    """
    Cut a smaller solvent box from a bigger one that will be used during the simulation.

    Given a (in the best case) larger solvent box, a smaller one will be cut.
    Having only as many solvent molecules as absolutely necessary reduces the
    calculation time. This is possible since the potential energy for a group
    of atoms may be calculated with lammps. The center of the box will be at
    the origin.

    > lmpdat            str; lammps data file of the system to cut from
    > box               Box; size of the box to cut out
    > lmpdat_out       str; name of new lammps data file with cut coordinates
    """
    md_sys = _read_system(lmpdat, dcd)

    # generate planes that describe the box
    cart_box = copy.deepcopy(box)
    cart_box.box_lmp2cart()

    # define planes which form the box to cut (plane vectors always need the same
    # origin or they will be shifted)
    plane_ab = agv.get_plane(cart_box.crt_a, cart_box.crt_b)
    plane_ca = agv.get_plane(cart_box.crt_c, cart_box.crt_a)
    plane_bc = agv.get_plane(cart_box.crt_b, cart_box.crt_c)

    # opposite planes to ab, ca and bc
    plane_ab_c = [-1 * i for i in agv.get_plane(cart_box.crt_a, cart_box.crt_b, cart_box.crt_c)]
    plane_ca_b = [-1 * i for i in agv.get_plane(cart_box.crt_c, cart_box.crt_a, cart_box.crt_b)]
    plane_bc_a = [-1 * i for i in agv.get_plane(cart_box.crt_b, cart_box.crt_c, cart_box.crt_a)]

    # cut plane
    md_sys.cut_shape(-1, True, plane_ab, plane_ab_c, plane_ca, plane_ca_b, plane_bc, plane_bc_a)

    # enlarge box a little so atoms at the edge will not clash due to pbc
    cart_box.box_cart2lat()
    cart_box.ltc_a += 4
    cart_box.ltc_b += 4
    cart_box.ltc_c += 4
    cart_box.box_lat2lmp()

    # def new box vectors by given box angles
    md_sys.ts_boxes[0] = cart_box

    md_sys.ts_boxes[0].lmp_xy = None
    md_sys.ts_boxes[0].lmp_xz = None
    md_sys.ts_boxes[0].lmp_yz = None

    md_sys.mols_to_grps()
    md_sys.change_indices(incr=1, mode="increase")
    md_sys.write_lmpdat(lmpdat_out, cgcmm=True)


# System Preparation ==========================================================#

def _rot_matrix():
    """
    Create an arbitrary rotational matrix.

    Returns
    -------
    m_arb_xyz : array of arrays
        The arbitrary rotational matrix

    """
    # rotation matrix, consisting of 3 arbitrary rotations A, B and C
    m_arb_rot_x = agm.arb_rot_matrix([1, 0, 0])  # A
    m_arb_rot_y = agm.arb_rot_matrix([0, 1, 0])  # B
    m_arb_rot_z = agm.arb_rot_matrix([0, 0, 1])  # C
    # multiply all three matrices; ABC=A(BC)=(AB)C
    m_arb_xy = np.matmul(m_arb_rot_x, m_arb_rot_y)  # AB
    m_arb_xyz = np.matmul(m_arb_xy, m_arb_rot_z)  # ABC
    return m_arb_xyz


def _rotate_sys(md_sys):
    """
    Rotate system arbitrarily.
    """
    m_arb_xyz = _rot_matrix()
    sys_add_all_atms = range(len(md_sys.atoms))
    md_sys.mm_atm_coords(-1, m_arb_xyz, False, *sys_add_all_atms)


def _shift_sys(md_sys, radius, radius_buffer=1):
    """
    """
    # enlarge radius so atoms do not clash
    radius += radius_buffer
    rn_pos = agm.points_on_sphere(npoints=1, ndim=3, radius=radius)[0]
    mx_trans = cgt.translation_matrix(rn_pos)
    sys_add_all_atms = range(len(md_sys.atoms))
    md_sys.mm_atm_coords(-1, mx_trans, False, *sys_add_all_atms)


def _create_new_box(md_sys):
    """
    """
    # add new orthogonal box
    md_sys.ts_boxes = []  # delete all previous unneeded boxes
    box_diameter = md_sys.get_system_radius(-1) + 30
    pi_2 = math.pi / 2
    new_box = mdb.Box(boxtype="lattice", ltc_a=box_diameter, ltc_b=box_diameter,
                      ltc_c=box_diameter, ltc_alpha=pi_2, ltc_beta=pi_2,
                      ltc_gamma=pi_2)
    md_sys.ts_boxes.append(new_box)
    md_sys.ts_boxes[0].box_lat2lmp(triclinic=False)


def _check_sys(md_sys, atm_idx_max):
    """
    """
    # check interatomic distances
    md_sys.create_linked_cells(-1)
    close_atms = md_sys.chk_atm_dist(frame_id=-1, min_dist=1.2)
    close_atms = sorted(close_atms)

    # check if new atoms were placed too near to the main system by comparing
    # the atom indices (atom indices lower than those from the original main
    # system are not counted)
    try:
        return close_atms[-1] > atm_idx_max
    except IndexError:
        return True


def sysprep(lmpdat_out, lmpdat_main, lmpdat_add, dcd_add=None,
            frame_idx=-1):
    """
    Prepare the system for the next docking step.

    Since the kawska-zahn approach adds an agglomerate of molecules each stage,
    the system must be prepared for lammps beforehand. Therefor a sphere with
    radius r1 is created around the main system as well as a sphere with radius
    r2 around the agglomerate to add. Both systems are combined, checked for clashes
    and eventually a new lammps data file is written which can further be
    utilized. If

    Parameters
    ----------
    lmpdat_out : str
        Name of merged lammps data file to write

    lmpdat_main : str
        Name of the lammps data file to add another system to

    lmpdat_add : str
        Name of the lammps data file which is added to the main molecular
        system

    Returns
    -------
    success : bool
        True if successful, False otherwise.

    Writes a new lammps data file called sysprep_out_index

    """
    # read and transpose the main sys to the origin
    main_sys = _read_system(lmpdat_main)
    main_sys.transpose_by_cog(-1, [0, 0, 0], copy=False)
    _natoms = len(main_sys.atoms)

    # read and transpose the add sys to the origin
    add_sys = _read_system(lmpdat_add, dcd_add, frame_idx)
    add_sys.transpose_by_cog(-1, [0, 0, 0], copy=False)

    # rotate add sys
    _rotate_sys(add_sys)

    # shift add sys to sphere around main sys
    main_sys_radius = main_sys.get_system_radius(-1)
    add_sys_radius = add_sys.get_system_radius(-1)
    kwz_radius = main_sys_radius + add_sys_radius

    _shift_sys(add_sys, kwz_radius)

    # merge both systems
    main_sys.extend_universe(add_sys, mode="merge")
    _create_new_box(main_sys)

    # group atoms by bonds
    main_sys.fetch_molecules_by_bonds()
    main_sys.mols_to_grps()

    success = _check_sys(main_sys, _natoms)

    # write an output lammps data only if everything worked out
    if success is True:
        # write new data file
        main_sys.change_indices(incr=1, mode="increase")
        main_sys.write_lmpdat(lmpdat_out, frame_id=0, title="System ready" +
                              "for docking", cgcmm=True)
    return success


# Quenching ===================================================================#
def quench(lmpcuts, lmpdat_main):
    """
    """
    main_sys = _read_system(lmpdat_main)
    natoms_main_sys = len(main_sys.atoms)
    del main_sys

    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp, neigh=True)

    lmp.file(lmpcuts.settings_file)
    # change box type to not be periodic
    lmp.command("boundary f f f")
    lmpcuts.read_system(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    # intermediate restart files
    lmp.command("restart {0} {1} {1}".format(lmpcuts.logsteps * 10, lmpcuts.output_lmprst))

    # define the atoms that may move during the simulation
    lmp.command("group grp_add_sys id > {}".format(natoms_main_sys))
    lmp.command("group grp_main_sys id <= {}".format(natoms_main_sys))
    #lmp.command("fix freeze grp_main_sys setforce {0} {0} {0}".format(0.0))

    # pre-optimization
    lmp.command("min_style cg")
    lmp.command("min_modify dmax 0.5")
    lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

    # set an additional push to the added atoms that should be docked
    # towards the origin at (0/0/0)
    # (prevents losing atoms due to being localized outside the cutoff)
    if rank == 0:
        prep_sys = _read_system(lmpcuts.input_lmpdat)
        cog = agm.get_cog(prep_sys.ts_coords[-1][natoms_main_sys + 1:])
        cog /= np.linalg.norm(cog, axis=0)  # unit vector
        cog *= -1
    else:
        cog = None

    cog = comm.bcast(cog, 0)
    lmp.command(("fix force grp_add_sys addforce {0} {1} {2} every 10000").format(*cog))

    # barostatting, thermostatting
    lmp.command("fix integrator grp_add_sys nvt temp {0} {1} 0.1".format(lmpcuts.tstart, lmpcuts.tstop))
    lmp.command("run {}".format(lmpcuts.runsteps))
    # remove pushing force
    lmp.command("unfix force")

    # post-optimization
    lmp.command("min_style quickmin")
    lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

    # write restart file
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()


# Annealing ===================================================================#
def relax_box(lmpcuts):
    """
    """
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp, neigh=True)

    lmp.file(lmpcuts.settings_file)
    lmpcuts.read_system(lmp)
    lmpcuts.thermo(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    # minimize cut box if not done already
    if lmpcuts.input_lmprst is None or os.path.isfile(lmpcuts.input_lmprst):
        lmpcuts.minimize(lmp, style="cg")

    # barostatting, thermostatting
    lmpcuts.berendsen(lmp)

    #lmp.command("fix nose_hoover grp_add_sys nvt temp {0} {1} 0.1".format(lmpcuts.tstart, lmpcuts.tstop))

    lmp.command("run {}".format(lmpcuts.runsteps))
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()


def fix_indent_ids(radii, cogs, group_name, scale_start=1, scale_stop=12):
    """
    Create the strings that are used as IDs for the lammps fix indent command.

    Parameters
    ----------
    group_name : str
        Name for the group, e.g. molecule or atom
    radii : list of floats
        Radii of all molecules to create a sphere around
    cogs : list of floats
        Centers of geometry of all molecules, i.e. center of each sphere

    Returns
    -------
    indent_fixes : list of lists
        Contains all ready up fixes for the lammps indent fix

    """
    if rank == 0:
        print(radii)

    fix_str = "fix {0} all indent 1 sphere {c[0]} {c[1]} {c[2]} {1} side out"
    all_indent_fixes = []

    for scaling_factor in xrange(scale_start, scale_stop):
        indent_fixes = []
        # set fixes for void_indent
        for idx, (radius, cog) in enumerate(zip(radii, cogs)):
            fix_name = "{}_{}_fix_indent".format(group_name, idx)
            grown_radius = radius * 0.1 * scaling_factor
            cur_fix_str = fix_str.format(fix_name, grown_radius, c=cog)
            indent_fixes.append(cur_fix_str)
        all_indent_fixes.append(indent_fixes)

    return all_indent_fixes


def lmp_indent(lmp, indents, runsteps, keep_last_fixes=False):
    """
    Perform indentation runs to grow spheres around certain molecules.

    Parameters
    ----------
    lmp : lammps-object
    indents : list of lists of str
    runsteps : int
    keep_last_fixes : bool

    """
    # make sphere grow (i.e. create and delete fixes for fix indent)
    #pdb.set_trace()

    for indent_fixes in indents[:-1]:
        for indent_fix in indent_fixes:
            lmp.command(indent_fix)
        lmp.command("run {}".format(runsteps))

        for indent_fix in indent_fixes:
            lmp.command("unfix {}".format(indent_fix.split()[1]))

    # last run for sphere growth
    for indent_fix in indents[-1]:
        lmp.command(indent_fix)
    lmp.command("run {}".format(runsteps))

    # unfix last fixes
    if keep_last_fixes is False:
        for indent_fix in indents[-1]:
            lmp.command("unfix {}".format(indent_fix.split()[1]))


def check_imaginary_clashes(sys_both, sys_a, sys_b, dcd_b):
    """
    Check for imaginary clashes between atoms of two systems before combining them.

    Parameters
    ----------
    sys_both : Universe
    sys_a : Universe
    sys_b : Universe

    Returns
    -------
    close_atms_a : list of int
        Atom ids of atoms that need additional spheres because they clash with
        solvent atoms otherwise

    """
    #sys_both.reset_cells()
    # read latest solvent coordinates and boxes
    sys_b.import_dcd(dcd_b)
    sys_b.read_frames(frame=-2, to_frame=-1)
    sys_b.close_dcd()

    # concatenate solute and latest solvent coordinates
    sys_both.ts_coords.append(np.concatenate((sys_a.ts_coords[-1], sys_b.ts_coords[-1])))
    sys_both.ts_boxes = sys_b.ts_boxes
    sys_both.create_linked_cells(-1, rcut_a=2, rcut_b=2, rcut_c=2)
    close_atms_ab = sys_both.chk_atm_dist(-1, min_dist=1.0, exclude_same_molecule=True)
    # only get close atoms from system a
    close_atms_a = [i for i in close_atms_ab if i <= len(sys_a.atoms)]

    sys_b.reset_cells()
    sys_both.reset_cells()
    return close_atms_a


def create_voids(lmpcuts, lmpdat_solvate, dcd_solvate):
    """
    """
    # solvate with molecule radii and cogs
    solvate_sys = _read_system(lmpdat_solvate, dcd_solvate)
    radii_mol, cogs_mol = _molecules_radii(solvate_sys)
    indent_strs = fix_indent_ids(radii_mol, cogs_mol, "molecule", scale_start=10, scale_stop=15)

    # load lammps and lammps settings
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp)

    lmp.file(lmpcuts.settings_file)
    lmpcuts.read_system(lmp)
    lmpcuts.thermo(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    lmpcuts.berendsen(lmp)
    #lmp.command("unfix integrator")
    #lmp.command("fix integrator all nve/limit 1.0")
    lmp_indent(lmp, indent_strs, lmpcuts.runsteps, keep_last_fixes=True)

    # solution system
    solvent_sys = _read_system(lmpcuts.input_lmpdat)
    solution_sys = merge_sys(solvate_sys, solvent_sys)
    solution_sys.reset_cells()

    # check solvate atoms with too close contacts to solvent atoms
    close_atoms = check_imaginary_clashes(solution_sys, solvate_sys, solvent_sys, lmpcuts.output_dcd)
    cogs_atoms = [solvate_sys.ts_coords[-1][i] for i in close_atoms]
    radii_atoms = [mde.elements_mass_radii[round(solvate_sys.atm_types[solvate_sys.atoms[i].atm_key].weigh, 1)] for i in close_atoms]

    # gather close atoms and remember which were close
    all_close_atoms = []
    all_radii_atoms = []
    all_cogs_atoms = []

    all_close_atoms.extend(close_atoms)
    all_radii_atoms.extend(radii_atoms)
    all_cogs_atoms.extend(cogs_atoms)

    factor_start = 10
    factor_stop = factor_start + 1

    if close_atoms != []:
        # move solvent molecules away from close solvate atoms
        for _ in xrange(5):
            indent_strs = fix_indent_ids(all_radii_atoms, all_cogs_atoms, "atom", scale_start=factor_start, scale_stop=factor_stop)
            lmp_indent(lmp, indent_strs, lmpcuts.runsteps, keep_last_fixes=False)
            close_atoms = check_imaginary_clashes(solution_sys, solvate_sys, solvent_sys, lmpcuts.output_dcd)

            # add new close atoms to present ones or stop indenting
            if close_atoms == []:
                break

            # add further close atoms to present ones
            for atm_idx in close_atoms:
                print(close_atoms)
                if atm_idx not in all_close_atoms:
                    all_close_atoms.append(atm_idx)
                    all_radii_atoms.append(mde.elements_mass_radii[round(solvate_sys.atm_types[solvate_sys.atoms[i].atm_key].weigh, 1)])
                    all_cogs_atoms.append(solvate_sys.ts_coords[-1][atm_idx])

            # dynamically grow sphere around atoms
            factor_start += 1
            factor_stop = factor_start + 1

    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()

    del (solvate_sys, solvent_sys, solution_sys)

    if rank == 0:
        solvate_sys = _read_system(lmpdat_solvate, dcd_solvate)
        solvent_sys = _read_system(lmpcuts.input_lmpdat, lmpcuts.output_dcd)
        solution_sys = merge_sys(solvate_sys, solvent_sys)
        solution_sys.change_indices(incr=1, mode="increase")
        solution_sys.write_lmpdat(lmpcuts.output_lmpdat, cgcmm=True)

    return close_atoms == []


def anneal_1(lmpcuts, lmpdat_solvate, dcd_solvate, lmpdat_solvent=None, dcd_solvent=None):
    """
    """
    pass
