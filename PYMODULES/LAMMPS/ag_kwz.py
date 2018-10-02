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
initial_velocity = "velocity all create {} 8324798 rot yes dist gaussian"
_sqrt_3 = math.sqrt(3)


def lmp_thermostat(lmp, logsteps):
    """
    Log thermodynamic data.
    """
    thermargs = ["step", "temp", "press", "vol", "density",
                 "cella", "cellb", "cellc", "cellalpha", "cellbeta", "cellgamma",
                 "etotal", "pe", "evdwl", "ecoul", "ebond", "eangle", "edihed",
                 "eimp", "enthalpy"]
    lmp.command("thermo_style custom " + " ".join(thermargs))
    lmp.command("thermo_modify lost warn flush yes")
    #lmp.command("thermo_modify line multi format float %g")
    lmp.command("thermo {}".format(logsteps))


class LammpsSettings(object):
    """
    """
    def __init__(self, tstart=None, tstop=None, pmin=None, pmax=None, logsteps=None, runsteps=None,
                 pc_file=None, settings_file=None, input_lmpdat=None,
                 input_lmprst=None, output_lmprst=None, output_lmplog=None,
                 output_dcd=None, output_name=None, gpu=False):
        """
        """
        self.tstart = tstart
        self.tstop = tstop
        self.pmin = pmin
        self.pmax = pmax
        self.logsteps = logsteps
        self.runsteps = runsteps
        self.pc_file = pc_file
        self.settings_file = settings_file
        self.input_lmpdat = input_lmpdat
        self.input_lmprst = input_lmprst
        self.output_lmplog = output_lmplog
        self.output_dcd = output_dcd
        self.output_lmprst = output_lmprst
        self.output_name = output_name
        self.gpu = gpu


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
def quench(lmpdat_main, lmp_settings):
    """
    """
    main_sys = _read_system(lmpdat_main)
    natoms_main_sys = len(main_sys.atoms)
    del main_sys

    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmp_settings.output_lmplog))

    # load gpu package
    if lmp_settings.gpu:
        # neighbor list building on gpu (currently) leads to segmentation faults
        #lmp.command("package gpu 1 neigh yes")
        lmp.command("package gpu 1 neigh no")  # DEBUGGING
        lmp.command("suffix gpu")

    # load general settings (units, boundary, dimension, atom_style, etc. pp.)
    lmp.file(lmp_settings.settings_file)

    # change box type to not be periodic
    lmp.command("boundary f f f")

    # read data and create an initial velocity or read the restart file from
    # a previous run
    if os.path.isfile(lmp_settings.input_lmprst) is True:
        lmp.command("read_restart {}".format(lmp_settings.input_lmprst))
    else:
        lmp.command("read_data {}".format(lmp_settings.input_lmpdat))
        #lmp.command(initial_velocity.format(temps[0]))

    # read pair coefficients from file if provided
    if lmp_settings.pc_file is not None:
        lmp.file(lmp_settings.pc_file)

    # trajectory
    lmp.command("dump trajectory all dcd {} {}".format(lmp_settings.logsteps,
                                                       lmp_settings.output_dcd))
    lmp.command("dump_modify trajectory unwrap yes")

    # intermediate restart files
    lmp.command("restart {0} {1} {1}".format(lmp_settings.logsteps * 10,
                                             lmp_settings.output_lmprst))

    lmp_thermostat(lmp, lmp_settings.logsteps)
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
        prep_sys = _read_system(lmp_settings.input_lmpdat)
        cog = agm.get_cog(prep_sys.ts_coords[-1][natoms_main_sys + 1:])
        cog /= np.linalg.norm(cog, axis=0)  # unit vector
        cog *= -1
    else:
        cog = None

    cog = comm.bcast(cog, 0)
    lmp.command(("fix force grp_add_sys addforce {0} {1} {2} every 10000").format(*cog))

    # barostatting, thermostatting
    lmp.command("fix nose_hoover grp_add_sys nvt temp {0} {1} 0.1".format(lmp_settings.tstart,
                                                                          lmp_settings.tstop))
    lmp.command("run {}".format(lmp_settings.runsteps))
    # remove pushing force
    lmp.command("unfix force")

    # post-optimization
    #lmp.command("min_style fire")
    lmp.command("min_style quickmin")
    #lmp.command("min_style cg")
    #lmp.command("min_modify dmax 0.5")
    lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

    # write restart file
    lmp.command("undump trajectory")
    #lmp.command("unfix freeze")
    lmp.command("unfix nose_hoover")
    lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmp_settings.output_lmprst))
    lmp.command("clear")
    lmp.close()


# Annealing ===================================================================#
def relax_box(lmp_settings):
    """
    """
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmp_settings.output_lmplog))

    # load gpu package
    if lmp_settings.gpu:
        # neighbor list building on gpu (currently) leads to segmentation faults
        lmp.command("package gpu 1 neigh yes")
        #lmp.command("package gpu 1 neigh no")  # DEBUGGING
        lmp.command("suffix gpu")

    lmp.file(lmp_settings.settings_file)

    if os.path.isfile(lmp_settings.input_lmprst) is True:
        lmp.command("read_restart {}".format(lmp_settings.input_lmprst))
    else:
        lmp.command("read_data {}".format(lmp_settings.input_lmpdat))
        lmp.command(initial_velocity.format(lmp_settings.tstart))

    if lmp_settings.pc_file is not None:
        lmp.file(lmp_settings.pc_file)

    lmp_thermostat(lmp, lmp_settings.logsteps)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")

    # trajectory
    lmp.command("dump trajectory all dcd {} {}".format(lmp_settings.logsteps,
                                                       lmp_settings.output_dcd))
    lmp.command("dump_modify trajectory unwrap yes")

    # minimize cut box if not done already
    if os.path.isfile(lmp_settings.input_lmprst) is False:
        lmp.command("min_style cg")
        lmp.command("min_modify dmax 0.5")
        lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

    # barostatting, thermostatting
    lmp.command("fix integrator all nve")
    lmp.command("fix thermostat all temp/berendsen {} {} 0.5".format(lmp_settings.tstart,
                                                                     lmp_settings.tstop))
    lmp.command("fix barostat all press/berendsen iso {} {} 50".format(lmp_settings.pstart,
                                                                       lmp_settings.pstop))
    #lmp.command("fix nose_hoover grp_add_sys nvt temp {0} {1} 0.1".format(lmp_settings.tstart,
    #                                                                      lmp_settings.tstop))

    lmp.command("run {}".format(lmp_settings.runsteps))
    lmp.command("undump trajectory")
    lmp.command("unfix ic_prevention")
    lmp.command("unfix integrator")
    lmp.command("unfix thermostat")
    lmp.command("unfix barostat")
    lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmp_settings.output_lmprst))
    lmp.command("clear")
    lmp.close()


def create_settings(lmp_ptr, lmp_settings):
    """
    """
    lmp_ptr.command("log {} append".format(lmp_settings.output_lmplog))

    # load gpu package
    if lmp_settings.gpu:
        # neighbor list building on gpu (currently) leads to segmentation faults
        lmp_ptr.command("package gpu 1 neigh yes")
        #lmp.command("package gpu 1 neigh no")  # DEBUGGING
        lmp_ptr.command("suffix gpu")

    lmp_ptr.file(lmp_settings.settings_file)

    if os.path.isfile(lmp_settings.output_lmprst) is True:
        lmp_ptr.command("read_restart {}".format(lmp_settings.output_lmprst))
    else:
        lmp_ptr.command("read_restart {}".format(lmp_settings.input_lmprst))

    if lmp_settings.pc_file is not None:
        lmp_ptr.file(lmp_settings.pc_file)

    lmp_thermostat(lmp_ptr, lmp_settings.logsteps)
    lmp_ptr.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")

    # trajectory
    lmp_ptr.command("dump trajectory all dcd {} {}".format(lmp_settings.logsteps,
                                                           lmp_settings.output_dcd))
    lmp_ptr.command("dump_modify trajectory unwrap yes")

    # barostatting, thermostatting
    lmp_ptr.command("fix integrator all nve")
    lmp_ptr.command("fix thermostat all temp/berendsen {} {} 0.5".format(lmp_settings.tstart,
                                                                         lmp_settings.tstop))
    lmp_ptr.command("fix barostat all press/berendsen iso {} {} 50".format(lmp_settings.pstart,
                                                                           lmp_settings.pstop))
    #lmp_ptr.command("fix nose_hoover grp_add_sys nvt temp {0} {1} 0.1".format(lmp_settings.tstart,
    #                                                                      lmp_settings.tstop))


def set_spheres_1(lmp, nmols, radii, cogs, runsteps):
    """
    """
    fix_names = ["void_indent_molecule_{}".format(idx) for idx in xrange(nmols)]
    fix = "fix {0} all indent 1 sphere {c[0]} {c[1]} {c[2]} {1} side out"

    #=============================================#
    # initial spheres around all solvate molecules
    #=============================================#

    for run in xrange(1, 12):
        # set fixes for void_indent
        for fix_name, sphere, cog in zip(fix_names, radii, cogs):
            grown_sphere = (sphere * 0.1 * run)
            molecule_void_fix = fix.format(fix_name, grown_sphere, c=cog)
            lmp.command(molecule_void_fix)

        lmp.command("run {}".format(runsteps))

        # delete all fixes except for the last run
        if run < 11:
            for fix_name in fix_names:
                lmp.command("unfix {0}".format(fix_name))


def check_clashes(sys_both, sys_a, sys_b, dcd_b):
    """
    """
    sys_both.reset_cells()
    # read latest solvent coordinates and boxes
    sys_b.import_dcd(dcd_b)
    sys_b.read_frames(frame=-2, to_frame=-1)
    sys_b.close_dcd()

    # concatenate solute and latest solvent coordinates
    sys_both.ts_coords.append(np.concatenate((sys_a.ts_coords[-1], sys_b.ts_coords[-1])))
    sys_both.ts_boxes = sys_b.ts_boxes
    sys_both.create_linked_cells(-1, rcut_a=2, rcut_b=2, rcut_c=2)
    close_atms = sys_both.chk_atm_dist(-1, min_dist=2.0, exclude_same_molecule=True)
    sys_both.reset_cells()
    return close_atms


def merge(sys_a, sys_b, pair_coeffs=False):
    """
    """
    sys_both = copy.deepcopy(sys_a)
    sys_both.extend_universe(sys_b, u1_frame_id=-1, u2_frame_id=-1, mode="merge")
    #sys_both.ts_boxes = []
    #sys_both.ts_coords = []

    if pair_coeffs is True:
        sys_both.mix_pair_types(mode="ii", mix_style="arithmetic")

    sys_both.fetch_molecules_by_bonds()
    sys_both.mols_to_grps()
    return sys_both


def create_voids(lmpdat_solvate, dcd_solvate, lmpdat_solvent, lmp_settings):
    """
    """
    solvate_sys = _read_system(lmpdat_solvate, dcd_solvate)
    nmols_solvate = len(solvate_sys.molecules)
    natms_solvate = len(solvate_sys.atoms)
    solvate_info = _molecules_radii(solvate_sys)
    radii_sphere = solvate_info[0]
    cogs_solvate = solvate_info[1]

    solvent_sys = _read_system(lmpdat_solvent)
    solution_sys = merge(solvate_sys, solvent_sys)
    solution_sys.reset_cells()

    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    create_settings(lmp, lmp_settings)
    set_spheres_1(lmp, nmols_solvate, radii_sphere, cogs_solvate, lmp_settings.runsteps)


    #==============================================#
    # additional spheres around close solvate atoms
    #==============================================#

    close_atoms_solvate = []
    addnames_void_fix = []
    atoms_coords_void = []
    atoms_radii_void = []
    atoms_fixes_void = []

    # 50 attempts to make stuff grow large enough
    for _ in xrange(50):

        #================================================#
        # interatomic distances (solvate - solvent) check
        #================================================#

        # array to test if any atoms still have close contacts
        current_close_atoms_solvate = []

        if rank == 0:
            close_atoms_solution = check_clashes()

            for close_atom in close_atoms_solution:
                # append only atoms from current check that are not in close_atoms_solvate already
                if (close_atom < natoms_solvate_sys and close_atom not in close_atoms_solvate):
                    close_atoms_solvate.append(close_atom)

                # append all close atoms from current check
                if close_atom < natoms_solvate_sys:
                    current_close_atoms_solvate.append(close_atom)

            del (close_atom, close_atoms_solution)

            #==========================================#
            # create additional void fixes if necessary
            #==========================================#

            for atom_solvate in close_atoms_solvate:

                #=========#
                # atom fix
                #=========#

                solvate_atom_name_fix = "void_indent_atom_{}".format(atom_solvate)

                if solvate_atom_name_fix not in atoms_fixes_void:
                    atoms_fixes_void.append(solvate_atom_name_fix)

                    #=================#
                    # atom coordinates
                    #=================#

                    solvate_atom_coords = solvate_sys.ts_coords[-1][atom_solvate]
                    atoms_coords_void.append(solvate_atom_coords)

                    #===========#
                    # atom radii
                    #===========#

                    type_atom = solvate_sys.atoms[atom_solvate].atm_key
                    weigh_atom = round(solvate_sys.atm_types[type_atom].weigh, 1)

                    try:
                        solvate_atom_radius = mde.elements_mass_radii[weigh_atom]
                    except KeyError:
                        solvate_atom_radius = 2.0  # dummy radius if element was not found by mass

                    solvate_atom_radius += 1.5  # buffer
                    atoms_radii_void.append(solvate_atom_radius)

        #===============================================#
        # generate additional fixes for additional voids
        #===============================================#

        close_atoms_solvate = comm.bcast(close_atoms_solvate, 0)
        atoms_coords_void   = comm.bcast(atoms_coords_void, 0)
        atoms_radii_void    = comm.bcast(atoms_radii_void, 0)
        atoms_fixes_void    = comm.bcast(atoms_fixes_void, 0)
        current_close_atoms_solvate = comm.bcast(current_close_atoms_solvate, 0)
        void_pylmp_fixes = []

        # grow spheres (i.e. create additional fixes if necessary, keep old ones, so some MD-runs)
        if len(current_close_atoms_solvate) > 0:
            for atom_sphere_run in xrange(1, 12):
                for atom_void_fix_name, atom_void_sphere, atom_void_coords in zip(atoms_fixes_void,
                                                                                  atoms_radii_void,
                                                                                  atoms_coords_void):
                    growing_atom_sphere = atom_void_sphere * 0.1 * atom_sphere_run

                    # check if additional fix already exists
                    if rank == 0:
                        # pylammps instructions only on rank 0!
                        for lammps_fix in void_pylmp.fixes:
                            # skip all non-indent fixes
                            if lammps_fix["style"] != "indent":
                                continue

                            # do nothing if indentation already exists
                            #if lammps_fix["name"] == atom_void_fix_name:
                            #    break

                        else:  # fix not defined in lammps
                            atom_void_fix = ("fix {0} all indent 1 sphere {c[0]} {c[1]} {c[2]} {1} side out").format(
                                atom_void_fix_name,
                                growing_atom_sphere,
                                c=atom_void_coords)
                            void_pylmp_fixes.append(atom_void_fix)
                            del atom_void_fix

                        del (lammps_fix)

                void_pylmp_fixes = comm.bcast(void_pylmp_fixes, 0)
                atom_void_fix_name = comm.bcast(atom_void_fix_name, 0)

                for lammps_fix in void_pylmp_fixes:
                    void_lmp.command(lammps_fix)

                # delete all current fixes from the list
                void_pylmp_fixes = []
                void_lmp.command("run {}".format(void_steps))

                #DEBUGGING
                #print(void_pylmp.fixes)
                #time.sleep(2)

                if atom_sphere_run < 11:
                    void_lmp.command("unfix {}".format(atom_void_fix_name))

            del (atom_void_sphere, atom_void_coords,
                 atom_void_fix_name, atom_sphere_run,
                 void_pylmp_fixes)

        else:
            break

    del (close_atoms_solvate, current_close_atoms_solvate)

    #if rank == 0:
    #    print(void_pylmp.fixes)
    #    time.sleep(1)

    # write final restart file
    void_lmp.command("undump void_dcd")
    void_lmp.command("unfix ic_prevention")
    void_lmp.command("write_restart {}".format(void_solv_rst))
    void_lmp.command("clear")
    void_lmp.close()
    del (void_lmp, void_pylmp)

    lmp.command("run {}".format(lmp_settings.runsteps))
    lmp.command("undump trajectory")
    lmp.command("unfix ic_prevention")
    lmp.command("unfix integrator")
    lmp.command("unfix thermostat")
    lmp.command("unfix barostat")
    lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmp_settings.output_lmprst))
    lmp.command("clear")
    lmp.close()
