
import pdb
import datetime
import os
#import copy
import shutil as sl
import re
#import argparse
import math
#import time
import numpy as np
#from natsort import natsorted
#import itertools as it
import scipy.stats
from mpi4py import MPI
from lammps import lammps, PyLammps
import Transformations as cgt
import md_elements as mde
import md_box as mdb
import md_universe as mdu
#import ag_unify_md as agum
import ag_geometry as agm
import ag_lammps as aglmp
import ag_lmplog as agl
#import ag_vectalg as agv
import ag_statistics as ags
import ag_plotting as agplot
#import vmd

#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator


#==============================================================================#
# Helper functions
#==============================================================================#
def generate_timestamp():
    """[summary]

    [description]
    """
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y%m%d%H%M%S")
    return timestamp


def rename(src, dst):
    """
    Force renaming by deleting a folder with the same name as dst.
    """
    try:
        os.rename(src, dst)
    except OSError:
        sl.rmtree(dst)
        os.rename(src, dst)


def get_natms(lmpdat):
    """[summary]

    [description]

    Parameters
    ----------
    lmpdat : {[type]}
        [description]

    Returns
    -------
    [type]
        [description]
    """
    with open(lmpdat, "r") as f_in:
        line = f_in.readline()
        while line != '':
            line = f_in.readline()
            if "atoms" in line:
                return int(line.split()[0])


def get_remaining_cycles(total_cycles):
    """
    Scan current folder names for numbers to get the current cycle.

    Returns
    -------
    remaining_cycles : int
        remaining cycles

    """
    def get_next_cycle():
        """
        """
        def get_finished_cycles():
            """
            """
            fc = []
            folders_pwd = ["{}/{}".format(pwd, i) for i in os.listdir(pwd) if os.path.isdir(i)]

            # get last cycle from directory
            for folder in folders_pwd:
                # skip all folders where anything went wrong
                if "fail" in folder:
                    continue

                #print(folder)
                cycle = re.match(r'.*?([0-9]+)$', folder).group(1)
                cycle = int(cycle)

                # avoid duplicates
                if cycle not in fc:
                    fc.append(cycle)

            return fc

        finished_cycles = get_finished_cycles()

        try:
            current_cycle = max(finished_cycles)
        except ValueError:
            current_cycle = 0

        last_requench_file = pwd + "/requench_{}/".format(current_cycle) + "requench_{}.dcd".format(current_cycle)

        if len(finished_cycles) >= 1:

            if os.path.isfile(last_requench_file) is True:
                next_cycle = current_cycle + 1
            else:
                next_cycle = current_cycle

        else:
            next_cycle = 0

        return (current_cycle, next_cycle, last_requench_file)

    pwd = os.getcwd()
    current_cycle, next_cycle, last_requench_file = get_next_cycle()
    total_cycles = list(range(total_cycles))
    idx_next_cycle = total_cycles.index(next_cycle)
    remain_cycles = total_cycles[idx_next_cycle:]
    return (remain_cycles, last_requench_file)

    #===========================#
    # Molecule to add by pattern
    #===========================#
    #id_pattern = 0
    #num_patterns = len(args.pa) - 1  # indices start with 0 so one less
    #total_cycles = []
#
    #for cycle in range(args.cycles):
    #    if id_pattern > num_patterns:
    #        id_pattern = 0
    #    total_cycles.append((cycle, args.pa[id_pattern]))
    #    id_pattern += 1


def create_folder(folder):
    """
    Create folder or skip creation if it already exists
    """
    try:
        os.mkdir(folder)
    except OSError:
        print("***Info: Folder {} already exists!".format(folder))


def write_to_log(string, filename="kwz_log"):
    """
    Write string 's' to log file named 'kwz.log'.
    """
    with open("kwz_log", "a") as kwz_log:
        kwz_log.write(string)

################################################################################
# Ensembles
################################################################################


def md_simulation(lmpcuts, group, style, ensemble, keyword_min=None, keyword=None, unwrap_dcd=True):
    """
    Perform an md simulation.
    """
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp, neigh=True)

    lmp.file(lmpcuts.settings_file)

    # change dielectric
    if lmpcuts.dielectric is not None:
        lmp.command("dielectric {}".format(lmpcuts.dielectric))

    lmpcuts.load_system(lmp)

    lmpcuts.thermo(lmp)
    lmpcuts.dump(lmp, unwrap=unwrap_dcd)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    # distribute the available cores; option must be set after loading pair coeffs!
    lmp.command("comm_style tiled")
    lmp.command("balance 1.0 rcb")

    # pre-minimize system before running the actual simulation
    if lmpcuts.input_lmprst is None or os.path.isfile(lmpcuts.input_lmprst):
        #lmp.command("min_modify line quadratic")
        lmpcuts.minimize(lmp, style="cg", keyword=keyword_min)

    # barostatting, thermostatting
    lmp.command("fix ic_prevention all momentum {} linear 1 1 1 angular rescale".format(lmpcuts.momentum_steps))

    # set the group
    if group != "all":
        lmp.command(group)
        groupname = group.split()[1]
    else:
        groupname = "all"

    if style.lower() == "berendsen":
        lmpcuts.fix_berendsen(lmp, groupname, ensemble, keyword)
    elif style.lower() == "nose_hoover":
        lmpcuts.fix_hoover(lmp, groupname, ensemble, keyword)
    elif style.lower() == "langevin":
        lmpcuts.fix_langevin(lmp, groupname, ensemble, keyword)
    else:
        raise Warning("Style not (yet) implemented!")

    lmp.command("reset_timestep 0")

    try:
        lmp.command("run {}".format(lmpcuts.runsteps))
    except:
        # prevent deadlock by not nicely ending the whole program if more
        # than one rank is used
        if size > 1:
            MPI.COMM_WORLD.Abort()

    #lmp.command("run {}".format(lmpcuts.runsteps))
    #lmpcuts.minimize(lmp, style="cg", keyword=keyword_min)
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()


def nose_hoover_md(lmpcuts, group="all"):
    """
    """
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp, neigh=False)

    lmp.file(lmpcuts.settings_file)

    # change dielectric
    if lmpcuts.dielectric is not None:
        lmp.command("dielectric {}".format(lmpcuts.dielectric))

    lmpcuts.load_system(lmp)
    lmpcuts.thermo(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp, unwrap=True)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    # minimize cut box if not done already
    if lmpcuts.input_lmprst is None or os.path.isfile(lmpcuts.input_lmprst):
        lmpcuts.minimize(lmp, style="cg")

    lmp.command("reset_timestep 0")

    try:
        lmp.command("run {}".format(lmpcuts.runsteps))
    except:
        # prevent deadlock by not nicely ending the whole program if more
        # than one rank is used
        if size > 1:
            MPI.COMM_WORLD.Abort()

    lmpcuts.unfix_undump(pylmp, lmp)
    #lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()


#==============================================================================#
# System Preparation
#==============================================================================#

def _molecules_radii(lmp_sys):
    """
    Get center of geometry and sphere of each solvate molecule which is placed inside the solvent box.

    Parameters
    ----------
    lmp_sys : Unification

    Returns
    -------
    radii : list of floats
        radius of each molecule

    cogs : list of lists
        coordinates of the center of geometry of each molecule

    """
    radii = []
    cogs = []

    for molecule in lmp_sys.molecules:
        sphere, cog = lmp_sys.get_mol_radius(-1, *molecule)
        radii.append(sphere)
        cogs.append(cog)

    return (radii, cogs)


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
    sys_add_all_atms = list(range(len(md_sys.atoms)))
    md_sys.mm_atm_coords(-1, m_arb_xyz, False, *sys_add_all_atms)


def _shift_sys(md_sys, radius, radius_buffer=1):
    """
    """
    # enlarge radius so atoms do not clash
    radius += radius_buffer
    rn_pos = agm.points_on_sphere(npoints=1, ndim=3, radius=radius)[0]
    mx_trans = cgt.translation_matrix(rn_pos)
    sys_add_all_atms = list(range(len(md_sys.atoms)))
    md_sys.mm_atm_coords(-1, mx_trans, False, *sys_add_all_atms)


def _create_new_box(md_sys):
    """
    """
    # add new orthogonal box
    md_sys.ts_boxes = []  # delete all previous unneeded boxes
    #box_diameter = md_sys.get_system_radius(-1) + 100
    # diameter with an additional size of 20 should suffice since it
    # is quite expensive for simulation runs with solvent
    #box_diameter = md_sys.get_system_radius(-1) + 50
    box_diameter = md_sys.get_system_radius(-1) + 200  # very large box which will be altered later anyways
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


def sysprep(lmpdat_out, lmpdat_main, lmpdat_add, dcd_main=None, dcd_add=None, frame_idx_main=-1, frame_idx_add=-1, main_radius_scale=1.5):
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

    dcd_main : str (optional, default: None)
        dcd file to load on top of lmpdat_main

    dcd_add : str (optional, default: None)
        dcd file to load on top of lmpdat_add

    frame_idx_main : int
        frame index of frame to add from dcd_main to lmpdat_main

    frame_idx_add : int
        frame index of frame to add from dcd_add to lmpdat_add

    main_radius_scale : int or float
        Radius to scale the main systems radius by

    Returns
    -------
    success : bool
        True if successful, False otherwise.

    Writes a new lammps data file called sysprep_out_index

    """
    # read and transpose the main sys to the origin
    main_sys = aglmp.read_lmpdat(lmpdat_main, dcd_main, frame_idx_main)
    main_sys.transpose_by_cog(-1, [0, 0, 0], copy=False)
    _natoms = len(main_sys.atoms)

    # read and transpose the add sys to the origin
    add_sys = aglmp.read_lmpdat(lmpdat_add, dcd_add, frame_idx_add)
    add_sys.transpose_by_cog(-1, [0, 0, 0], copy=False)

    # rotate add sys
    _rotate_sys(add_sys)

    # shift add sys to sphere around main sys
    main_sys_radius = main_sys.get_system_radius(-1)
    add_sys_radius = add_sys.get_system_radius(-1)

    # old procedure
    #kwz_radius = main_sys_radius + add_sys_radius

    # new procedure: scale radius by main_sys size
    kwz_radius = main_sys_radius * radius_scale + add_sys_radius

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


################################################################################
# Quenching
################################################################################

def quench(lmpcuts, lmpdat_main, runs=20, split=None):
    """
    """
    # check how many cores are used, leave the list empty if it is only one
    try:
        other_ranks = list(range(lmpcuts.ncores))[1:]
    except ValueError:
        other_ranks = []

    def _check_success():
        """
        Check aggregation state and close lammps instance properly.
        """
        if rank == 0:
            quench_sys = aglmp.read_lmpdat(lmpcuts.input_lmpdat, dcd=lmpcuts.output_dcd)
            succeeded = quench_sys.check_aggregate()

            # send data to the other ranks
            for other_rank in other_ranks:
                comm.send(succeeded, dest=other_rank)

        else:
            succeeded = comm.recv(source=0)

        #print(succeeded)
        #succeeded = comm.bcast(succeeded, 0)

        # stop trying if it was
        if succeeded is True:
            # write restart file
            lmpcuts.unfix_undump(pylmp, lmp)
            lmp.command("reset_timestep 0")
            lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
            lmp.command("clear")
            lmp.close()

        return succeeded

    def _run(steps):
        """
        Helper function for running and catching an exception when anything
        goes wrong.
        """
        try:
            lmp.command("run {}".format(steps))
        except:
            # prevent deadlock by not nicely ending the whole program if more
            # than one rank is used
            if size > 1:
                MPI.COMM_WORLD.Abort()

    natoms_main_sys = get_natms(lmpdat_main)

    lmp = lammps(comm=split)
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp, neigh=False)

    lmp.file(lmpcuts.settings_file)

    # change dielectric
    if lmpcuts.dielectric is not None:
        lmp.command("dielectric {}".format(lmpcuts.dielectric))

    lmpcuts.load_system(lmp)

    #lmp.command("velocity all create {} {} mom yes rot yes dist gaussian".format(lmpcuts.tstart, np.random.randint(29847587)))
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp, unwrap=True)
    lmpcuts.thermo(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    # distribute the available cores if we have lots of vacuum in our box
    lmp.command("comm_style tiled")
    lmp.command("balance 1.0 rcb")

    # define the atoms that may move during the simulation
    lmp.command("group grp_add_sys id > {}".format(natoms_main_sys))
    #lmp.command("group grp_main_sys id <= {}".format(natoms_main_sys))
    #lmp.command("fix freeze grp_main_sys setforce {0} {0} {0}".format(0.0))

    # pre-optimization
    lmp.command("min_style cg")
    lmp.command("min_modify dmax 0.5")
    lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

    # set an additional push to the added atoms that should be docked
    # towards the origin at (0/0/0)
    # (prevents losing atoms due to being localized outside the cutoff)
    if rank == 0:
        prep_sys = aglmp.read_lmpdat(lmpcuts.input_lmpdat)
        cog = agm.get_cog(prep_sys.ts_coords[-1][natoms_main_sys + 1:])
        cog /= np.linalg.norm(cog, axis=0)  # unit vector
        # make vector show towards the center (0/0/0)
        cog_force = cog * -1

        for other_rank in other_ranks:
            comm.send(cog_force, dest=other_rank)

    else:
        cog_force = comm.recv(source=0)

    #cog_force = comm.bcast(cog_force, 0)
    # barostatting, thermostatting only for atoms that will be docked
    lmp.command("fix integrator grp_add_sys nvt temp {0} {1} 0.1".format(lmpcuts.tstart, lmpcuts.tstop))
    quench_success = False

    # runs attempts to dock the molecule
    for _ in range(runs):
        # minimize and check if that is enough for docking
        lmp.command("min_style quickmin")
        lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

        lmp.command("min_style cg")
        lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

        # check aggregate, i.e. docking was a success
        quench_success = _check_success()

        if quench_success is True:
            break

        del quench_success
        # perform a little molecular dynamics simulation with
        # half of the lmpcuts.runsteps
        _run(int(lmpcuts.runsteps))

        # check aggregate, i.e. docking was a success
        quench_success = _check_success()

        if quench_success is True:
            break

        del quench_success
        # give 'to-be-docked' molecules a little push if minimization was not
        # sufficient for aggregation
        addforce_cmd = "fix push grp_add_sys addforce {c[0]} {c[1]} {c[2]} every 1"
        addforce_cmd = addforce_cmd.format(c=cog_force)
        lmp.command(addforce_cmd)
        #import time
        #time.sleep(30)

        # set force does not work
        #lmp.command("fix push grp_add_sys setforce {0} {0} {0}".format(*cog))
        _run(1)

        # remove pushing force
        lmp.command("unfix push")

        # run the 2nd half of the simulation with push 'grp_add_sys'
        _run(int(lmpcuts.runsteps * 0.5))

        # stop to-be-docked molecules from moving
        #lmp.command("velocity grp_add_sys set 0.0 0.0 0.0")
        #lmp.command("fix freeze grp_add_sys setforce {0} {0} {0}".format(0.0))

        _run(int(lmpcuts.runsteps * 0.25))

    return quench_success


################################################################################
# Annealing
################################################################################

def _fix_indent_ids(radii, cogs, group_name, scale_start=1, scale_stop=12):
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
    fix_str = "fix {0} all indent 1 sphere {c[0]} {c[1]} {c[2]} {1} side out"
    all_indent_fixes = []

    for scaling_factor in range(scale_start, scale_stop):
        indent_fixes = []
        # set fixes for void_indent
        for idx, (radius, cog) in enumerate(zip(radii, cogs)):
            fix_name = "{}_{}_fix_indent".format(group_name, idx)
            grown_radius = radius * 0.1 * scaling_factor
            cur_fix_str = fix_str.format(fix_name, grown_radius, c=cog)
            indent_fixes.append(cur_fix_str)
        all_indent_fixes.append(indent_fixes)

    return all_indent_fixes


def _lmp_indent(lmp, indents, runsteps, keep_last_fixes=False):
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
    if not keep_last_fixes:
        for indent_fix in indents[-1]:
            lmp.command("unfix {}".format(indent_fix.split()[1]))


def _check_clashes(sys_both, sys_a, sys_b, dcd_b=None, unwrap=False):
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
    if dcd_b:
        sys_b.import_dcd(dcd_b)
        sys_b.read_frames(frame=-2, to_frame=-1)
        sys_b.close_dcd()

    # caveat: unwrap only if totally necessary; if center(box) is (0,0,0) and
    # cog of the molecule is also (0,0,0) then clashes can safely be measured
    if unwrap:
        sys_b.unwrap_cell(frame_id=-1)

    # testing #################################################################
    #sys_b.change_indices()
    #sys_b.write_lmpdat("unwrapped.lmpdat", frame_id=-1, cgcmm=True)
    #exit()
    # #########################################################################

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


def create_voids(lmpcuts, lmpdat_solvate, dcd_solvate=None, dcd_solvent=None):
    """
    """
    # solvate with molecule radii and cogs
    solvate_sys = aglmp.read_lmpdat(lmpdat_solvate, dcd_solvate, frame_idx_start=-2, frame_idx_stop=-1)
    radii_mol, cogs_mol = _molecules_radii(solvate_sys)
    indent_strs = _fix_indent_ids(radii_mol, cogs_mol, "molecule", scale_start=10, scale_stop=15)

    # load solution system (last frame only)
    solvent_sys = aglmp.read_lmpdat(lmpcuts.input_lmpdat, dcd_solvent, frame_idx_start=-2, frame_idx_stop=-1)
    solution_sys = mdu.merge_systems([solvate_sys, solvent_sys])
    solution_sys.reset_cells()

    # check solvate atoms with too close contacts to solvent atoms
    close_atoms = _check_clashes(solution_sys, solvate_sys, solvent_sys)
    cogs_atoms = [solvate_sys.ts_coords[-1][i] for i in close_atoms]
    radii_atoms = [mde.elements_mass_radii[round(solvate_sys.atm_types[solvate_sys.atoms[i].atm_key].weigh, 1)] for i in close_atoms]

    # gather close atoms and remember which were close
    all_close_atoms = []
    all_radii_atoms = []
    all_cogs_atoms = []

    all_close_atoms.extend(close_atoms)
    all_radii_atoms.extend(radii_atoms)
    all_cogs_atoms.extend(cogs_atoms)

    # load lammps and lammps settings
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {}".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp)

    lmp.file(lmpcuts.settings_file)

    # change dielectric
    if lmpcuts.dielectric is not None:
        lmp.command("dielectric {}".format(lmpcuts.dielectric))

    lmpcuts.load_system(lmp)
    lmpcuts.thermo(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp, unwrap=True)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    #lmpcuts.fix_berendsen(lmp, group="all", ensemble="nve", keyword="iso")
    lmpcuts.fix_berendsen(lmp, group="all", ensemble="npt", keyword="iso", integrator="nve/limit 0.2")
    #lmp.command("unfix integrator")
    #lmp.command("fix limit_movement all nve/limit 0.05")
    _lmp_indent(lmp, indent_strs, lmpcuts.runsteps, keep_last_fixes=True)

    factor_start = 10
    factor_stop = factor_start + 1

    if close_atoms != []:
        # move solvent molecules away from close solvate atoms
        for _ in range(5):
            indent_strs = _fix_indent_ids(all_radii_atoms, all_cogs_atoms, "atom", scale_start=factor_start, scale_stop=factor_stop)
            _lmp_indent(lmp, indent_strs, lmpcuts.runsteps, keep_last_fixes=False)
            close_atoms = _check_clashes(solution_sys, solvate_sys, solvent_sys, lmpcuts.output_dcd)

            # add new close atoms to present ones or stop indenting
            if close_atoms == []:
                break

            # add further close atoms to present ones
            for atm_idx in close_atoms:
                #print(close_atoms)
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

    return close_atoms == []


def requench(lmpcuts, minstyle="cg"):
    """[summary]

    [description]

    Parameters
    ----------
    lmpcuts : {[type]}
        [description]
    minstyle : {str}, optional
        [description] (the default is "cg", which [default_description])

    Returns
    -------
    [type]
        [description]
    """
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    lmp.file(lmpcuts.settings_file)

    # change dielectric
    if lmpcuts.dielectric is not None:
        lmp.command("dielectric {}".format(lmpcuts.dielectric))

    # read restart file from previous run
    lmpcuts.load_system(lmp)
    lmpcuts.thermo(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmp.command("fix integrator all nvt temp {0} {1} 0.1".format(lmpcuts.tstart, lmpcuts.tstop))
    lmpcuts.dump(lmp, unwrap=True)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    lmp.command("reset_timestep 0")
    lmp.command("run {}".format(lmpcuts.runsteps))
    lmpcuts.minimize(lmp, style=minstyle)
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()

    # check if aggregate is intact after requenching
    md_sys = aglmp.read_lmpdat(lmpcuts.input_lmpdat, lmpcuts.output_dcd)
    aggregate_ok = md_sys.check_aggregate(frame_id=-1)
    return aggregate_ok


def _append_data(data, lmplog, fstart=1, thermo="PotEng"):
    """
    fstart : int
        frame to start recording from
    """
    cur_log = agl.LmpLog()
    cur_log.read_lmplog(lmplog)

    # potential energy of the whole system (solvent included if part of system)
    cur_data = cur_log.data[-1][thermo][fstart:]
    data.extend(cur_data)


#def vmd_rmsd_and_cluster(lmpdat_solution, lmpdat_solvate, dcd_files, atm_idxs, xyz_out, percentage_to_check=80):
#    """
#    Align solvate molecules and calculate the rmsd.
#
#    The first frame serves as the reference frame for the alignment and the rmsd
#    calculation. The vmd modules molecule and atomsel have to be loaded for
#    this function to work (and the coordinates had also to be read already).
#
#    Parameters
#    ----------
#    lmpdat_solution : str
#        name of lammps data file with topology to load
#    dcd_files : list of str
#        file names of the dcd files to load into vmd
#
#    Returns
#    -------
#    rmsds : list of floats
#        rmsd values for each aligned selection
#    cluster0 : list of ints
#        indices of frames from the most populated cluster
#
#    """
#    import vmd
#    import molecule
#    import atomsel
#
#    # get number of solvate atoms
#    natms_solvate = get_natms(lmpdat_solvate)
#
#    # get number of frames to check
#    percentage_to_check /= 100.0
#    # read topology and coordinates
#    molecule.load("lmpdat", lmpdat_solution)
#
#    for dcd_file in dcd_files:
#        molecule.read(0, "dcd", dcd_file, waitfor=-1)
#
#    # rmsd with alignment for the solvent atoms
#    rmsds = []
#    selection = "index {}".format(" ".join(map(str, atm_idxs)))
#    sel1 = atomsel.atomsel(selection, frame=0)
#
#    frame_first = int(1 - percentage_to_check * molecule.numframes(0))
#    frame_last = molecule.numframes(0)
#
#    for frame_idx in range(frame_first, frame_last):
#        sel2 = atomsel.atomsel(selection, frame=frame_idx)
#        matrix_algin = sel2.fit(sel1)
#        sel2.move(matrix_algin)
#        rmsd = sel2.rmsd(sel1)
#        rmsds.append(rmsd)
#
#    # clustering (currently only possible using the tcl interface?)
#    tcl_selection = 'set tcl_sel1 [atomselect {} "index 0 to {} frame {}"]'
#    tcl_sel1 = tcl_selection.format(0, natms_solvate - 1, 0)
#    vmd.evaltcl(tcl_sel1)
#    cluster = "measure cluster $tcl_sel1 num 5 distfunc rmsd cutoff 1.0 first {} last {} step 1 selupdate True"
#    clusters = vmd.evaltcl(cluster.format(frame_first, frame_last))
#    clusters = re.findall(r"\{([^[\}]*)\}", clusters)
#    cluster_0 = [int(i) for i in clusters[0].split()]
#    molecule.write(0, "xyz", xyz_out, beg=cluster_0, end=cluster_0, selection=sel1)
#    vmd.VMDexit("closing vmd")
#    return (rmsds, cluster_0)


################################################################################
# Annealing
################################################################################

def _anneal(lmpcuts, pe_atm_idxs, ensemble, group="all", keyword="iso"):
    """Helper function for the annealing step of the kawska-zahn approach.

    This function calls lammps to execute the simulation for the annealing
    part of the kawska-zahn approach.

    Parameters
    ----------
    lmpcuts : {[type]}
        [description]
    pe_atm_idxs : {list of int}
        List of atom indices for the solvate atoms. Needed to compute the potential
        energy of the solvate atoms only.
    ensemble : {[type]}
        [description]
    group : {str}, optional
        [description] (the default is "all", which [default_description])
    keyword : {str}, optional
        [description] (the default is "iso", which [default_description])
    """
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    # caveat: only possible with neigh no if calculating pe/atom on the graphics card
    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp, neigh=False)

    lmp.file(lmpcuts.settings_file)

    # change dielectric
    if lmpcuts.dielectric is not None:
        lmp.command("dielectric {}".format(lmpcuts.dielectric))

    lmpcuts.load_system(lmp)

    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp, unwrap=True)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    # distribute the available cores if we have lots of vacuum in our box
    lmp.command("comm_style tiled")
    lmp.command("balance 1.0 rcb")

    # compute potential energy of the solvate (pair, bond, ...)
    # in order to compute the pe of the solvent, the solvent atom-ids must be
    # known
    pe_atm_ids = [i + 1 for i in pe_atm_idxs]
    atm_ids_str = " ".join(map(str, pe_atm_ids))
    lmp.command("group resname_atoms id {}".format(atm_ids_str))
    lmp.command("compute pe_per_atom_solvate resname_atoms pe/atom")
    lmp.command("compute pe_solvate_complete resname_atoms reduce sum c_pe_per_atom_solvate")

    if "c_pe_solvate_complete" not in lmpcuts.thermargs:
        lmpcuts.thermargs.append("c_pe_solvate_complete")

    lmpcuts.thermo(lmp, hb_group="resname_atoms")
    lmpcuts.fix_hoover(lmp, group, ensemble, keyword)

    # reset the timestep in order to have a consistent number of steps in the
    # log- and dcd-file (steps may differ after e.g. a minimization)
    lmp.command("reset_timestep 0")

    lmp.command("run {}".format(lmpcuts.runsteps))
    lmpcuts.unfix_undump(pylmp, lmp)
    # not sure if resetting the timestep is necessary
    #

    # write restart file only if equilibrated (so definitely not here!)
    lmp.command("write_restart {}".format(lmpcuts.inter_lmprst))
    lmp.command("clear")
    lmp.close()


def _test_anneal_equil(data, output=None, xlabel=None):
    """
    Check normality distribution of the underlying data.

    If the qq-plot fails, the distribution is not normal, if the other tests
    fail, normality can still be accepted if the qq plot states so.
    The plots of the data (histogram, QQ-plot) will be saved in order to quickly
    check the distribution manually.

    Parameters
    ----------
    data : list of floats
        Values to check (1-dimensional array).
    output : str
        Prefix of all output files.
    xlabel : str
        Name of the xlabel for the histogram plot.

    Returns
    -------
    equilibrated : bool
        Result if the values show normal distribution. Normal distribution is
        only accepted if the linear regression of the QQ-plot has a correlation-
        coefficient larger than 0.998, the curve is not too skewed.

    """
    # save the histogram plot
    agplot.gnuplot_gaussfit_plot(data, xlabel=xlabel, output=output + "_histogram")

    qq_normal = ags.qq_test(data, rsquare_thrsh=0.998, output=output + "_qq_plot", save_plot=True)
    skew_normal = ags.test_gauss_shape("skewness", data)

    try:
        kurtosis_normal = ags.test_gauss_shape("kurtosis", data)
    except Warning:
        kurtosis_normal = False

    # does not work properly with many
    equilibrated = qq_normal and (skew_normal or kurtosis_normal)
    print("QQ-Normal: {}, Kurtosis: {}, Skew: {}".format(qq_normal, kurtosis_normal, skew_normal))

    return equilibrated


def anneal_productive(lmpcuts, atm_idxs_solvate, percentage_to_check, ensemble, group="all", keyword=None, attempts=500, output=None):
    """
    Carry out a productive run for the annealing step of the Kawska-Zahn approach.

    Check regularly if the simulation is equilibrated (i.e. PotEng of the system
    is normally distributed). If that criteria is fulfilled, stop and check if
    the aggregate of the last frame is still intact (what is happening in between
    is not of our concern and will be checked later when the best conformation
    of the complex is chosen).

    Parameters
    ----------
    lmpcuts : ag_lammps_sim.LmpSim
        The simulation settings.
    atm_idxs_solvate : list or tuple
        Indices of the solvate atoms.
    percentage_to_check : int or float, float will be converted to int
        The percentage of all frames (starting from the last frame to the first)
        to check the equilibration state from.
    ensemble : str
        The ensemble to run the simulation with
    group : str, optional
        The group of atoms that will be integrated. Uses lammps syntax for
        grouping. The default value is 'all'
    keyword : str, optional
        The keyword for the box relaxation. Allowed keywords are the same as
        in the lammps manual for barostatting: iso, aniso, tri. Default is 'None',
        which means no barostatting is done (i.e. nvt).
    output : str, optional
        Prefix for all files which were/are written by the run.

    Returns
    -------
    aggregate_ok : boolean
        The status of the aggregate. 'True' if aggregate did not dissolve, during
        the run, 'False' if it did.
    dcd_files : list
        A list of all DCD files that were written during the run. For each sub-
        run, an int starting from '0' is prepended to the filename until the
        whole simulation (formed by all sub-runs) shows convergence.
    log_files : list
        A list of all lammps log files that were written during the run. For each sub-
        run, an int starting from '0' is prepended to the filename until the
        whole simulation (formed by all sub-runs) shows convergence. Log-files
        are always in accordance to the current DCD file.

    """
    all_data = []
    solvate_sys_natoms = len(atm_idxs_solvate)
    dcd_files = []
    log_files = []

    # renaming the current found files independently of their original name automatically
    # sorts them in the right order
    for run_idx in range(attempts):
        # get just the base name of the files
        dcd_path, dcd_filename = os.path.split(lmpcuts.output_dcd)
        log_path, log_filename = os.path.split(lmpcuts.output_lmplog)

        # create a filename having 'run_idx' and check if this file already exists
        run_idx_pattern = re.compile(r'^[0-9]+')

        # check if filename starts with an integer
        if run_idx_pattern.match(dcd_filename) is not None:
            # split dcd-file by '_'
            split_dcd_filename = dcd_filename.split("_")
            dcd_filename = "{}_{}_{}".format(run_idx, split_dcd_filename[1], split_dcd_filename[2])
            del split_dcd_filename
        else:
            dcd_filename = str(run_idx) + "_" + dcd_filename

        if run_idx_pattern.match(log_filename) is not None:
            split_log_filename = log_filename.split("_")
            log_filename = "{}_{}_{}".format(run_idx, split_log_filename[1], split_log_filename[2])
            del split_log_filename
        else:
            log_filename = str(run_idx) + "_" + log_filename

        lmpcuts.output_dcd = "{}/{}".format(dcd_path, dcd_filename)
        lmpcuts.output_lmplog = "{}/{}".format(log_path, log_filename)

        # read previous log- and dcd-files if they exist already
        if os.path.isfile(lmpcuts.output_dcd) and os.path.isfile(lmpcuts.output_lmplog):
            dcd_files.append(lmpcuts.output_dcd)
            log_files.append(lmpcuts.output_lmplog)

            # read potential energy of all coordinates
            if rank == 0:

                # if this fails, then a lmplog was written but no simulation
                # was carried out
                try:
                    _append_data(all_data, lmpcuts.output_lmplog)
                except IndexError:
                    os.remove(lmpcuts.output_dcd)
                    os.remove(lmpcuts.output_lmplog)
                    dcd_files.pop(-1)
                    log_files.pop(-1)

            # skip the rest as often as there is no new file
            continue

        else:
            # carry out annealing run
            _anneal(lmpcuts, atm_idxs_solvate, ensemble, group, keyword)
            dcd_files.append(lmpcuts.output_dcd)
            log_files.append(lmpcuts.output_lmplog)

            # append data of last run
            if rank == 0:
                _append_data(all_data, lmpcuts.output_lmplog)

        # test given data so far (test applies to each newly carried out run)
        if rank == 0:
            # caveat:   each first frame of output_lmplog is omitted since
            #           it is not part of the dcd file
            # check last N % of all frames
            num_frames_to_check = int(percentage_to_check / 100 * len(all_data))

            # last X % of all frames (from end to start)
            normally_dstributed = _test_anneal_equil(all_data[-num_frames_to_check:], xlabel="Potential Energy / eV", output=output)

            # check if aggregate is still fine after the last run (only if we have a normal distribution)
            solution_sys = aglmp.read_lmpdat(lmpcuts.input_lmpdat, lmpcuts.output_dcd)
            solution_sys_atoms_idxs = list(range(len(solution_sys.atoms)))
            aggregate_ok = solution_sys.check_aggregate(frame_id=-1, excluded_atm_idxs=solution_sys_atoms_idxs[solvate_sys_natoms:])
            del (num_frames_to_check, solution_sys, solution_sys_atoms_idxs)
        else:
            aggregate_ok = False
            normally_dstributed = False

        aggregate_ok = comm.bcast(aggregate_ok)
        normally_dstributed = comm.bcast(normally_dstributed)

        # stop further runs if the aggregate is not ok or if the aggregate
        # is ok and the run is equilibrated (i.e. normal distribution of the
        # potential energy of the whole system)
        if (aggregate_ok is True and normally_dstributed is True) or aggregate_ok is False:
            break

        # the following prevents the 'else' condition (from for loop above) to execute
        #if run_idx == (attempts - 1):
        #    break

    # this only applies if all attempts already exist, but the
    # 'lmpcuts.output_lmprst' file was never written (e.g. due to aborted runs)
    else:
        aggregate_ok = False

        if rank == 0:
            num_frames_to_check = int(percentage_to_check / 100 * len(all_data))
            normally_dstributed = _test_anneal_equil(all_data[-num_frames_to_check:], xlabel="Potential Energy / eV", output=output)

            if normally_dstributed is True:
                solution_sys = aglmp.read_lmpdat(lmpcuts.input_lmpdat, lmpcuts.output_dcd)
                solution_sys_atoms_idxs = list(range(len(solution_sys.atoms)))
                aggregate_ok = solution_sys.check_aggregate(frame_id=-1, excluded_atm_idxs=solution_sys_atoms_idxs[solvate_sys_natoms:])
        else:
            normally_dstributed = False
            aggregate_ok = False

        aggregate_ok = comm.bcast(aggregate_ok, root=0)
        normally_dstributed = comm.bcast(normally_dstributed, root=0)

        #if aggregate_ok is True and normally_dstributed is True:
        #    sl.copy(lmpcuts.inter_lmprst, lmpcuts.output_lmprst)

    # if everything worked, name intermediate restart file as final output file
    if aggregate_ok is True and normally_dstributed is True:
        sl.copy(lmpcuts.inter_lmprst, lmpcuts.output_lmprst)

    return aggregate_ok and normally_dstributed


def find_best_frame(lmplogs, dcds, thermo="c_pe_solvate_complete", percentage_to_check=80):
    """
    Find the frame with the lowest energy / value for use in the requenching procedure.

    Parameters
    ----------
    lmplogs : tuple or list
        lammps-log files with thermo information.

    dcds : tuple or list
        dcd files with coordinates which correspond the given data of the
        log files (same number of frames is mandatory)

    Returns
    -------
    total_min_dcd : str
        name of dcd file with the lowest energy / value

    total_min_idx : int
        index of frame with the lowest energy / value

    total_min_val : int or float
        lowest value found over all lmplogs for thermo

    """
    #total_min_val = 1e20
    #total_min_idx = None
    #total_min_dcd = ""
    total_data = []

    # read the whole dataset, i.e. all lmplog files
    for lmplog in lmplogs:
        cur_log = agl.LmpLog()
        cur_log.read_lmplog(lmplog)

        # flatten list if more than one run was done
        # neglect the first frame since it is redundant in the log file
        cur_data = [i[thermo][1:] for i in cur_log.data]
        cur_data = [item for i in cur_data for item in i]
        total_data.append(cur_data)

    del cur_data

    # find minimum value in the last X % of the data
    total_data_flattened = [item for i in total_data for item in i]
    last_frames_to_check = int(percentage_to_check / 100 * len(total_data_flattened))
    data_to_check = total_data_flattened[-last_frames_to_check:]
    min_val = min(data_to_check)

    for cur_data, dcd in zip(total_data, dcds):
        if min_val in cur_data:
            return (dcd, cur_data.index(min_val), min_val)
    else:
        raise Exception("This should not have happened :/")

    #return (total_min_dcd, total_min_idx, total_min_val)


def write_requench_data(lmpdat_a, dcd_ab, index,
                        lmpdat_b=None,
                        output_lmpdat_a="output_name_a.lmpdat",
                        output_lmpdat_b="output_name_b.lmpdat"):
    """
    Write a new data file with the coordinates from a dcd file given by the index.

    Writes a new lammps data file, which has the same topology but different
    coordinates than the lmpdat that is given.

    Parameters
    ----------
    lmpdat_a : str
        lammps data file with the topology e.g. for the solvate

    lmpdat_b : None or str
        lammps data file with the topology e.g. for the solvent

    dcd_ab : str
        dcd file with frames to read from the solvate-solvent system

    index : int
        index to extract the frame from (i.e. index of energetically best frame)

    """
    sys_lmpdat_a = aglmp.read_lmpdat(lmpdat_a)
    sys_lmpdat_a_natoms = len(sys_lmpdat_a.atoms)

    sys_dcd_ab = aglmp.LmpStuff()
    sys_dcd_ab.import_dcd(dcd_ab)

    # read only the best frame which will be appended to the existing ones
    sys_dcd_ab.read_frames(frame=index, to_frame=index + 1)

    # apply box from dcd, since we are dealing with coordinates considering
    # that box (not doing this leads to errors during wrapping)
    sys_lmpdat_a.ts_boxes = sys_dcd_ab.ts_boxes

    # write only relevant coordinates for system a
    sys_lmpdat_a.ts_coords.append(sys_dcd_ab.ts_coords[-1][:sys_lmpdat_a_natoms])

    # wrap coordinates outside the box back inside (may happen when a molecule
    # flies of and reconnects with the aggregate)
    sys_lmpdat_a.wrap_cell(frame_id=-1, same_molecule=True)

    sys_lmpdat_a.change_indices()
    sys_lmpdat_a.write_lmpdat(output_lmpdat_a, -1, title="Best frame of {} with index {}".format(os.path.basename(dcd_ab), index), cgcmm=True)

    # write only relevant coordinates for system b
    if lmpdat_b is not None:
        sys_lmpdat_b = aglmp.read_lmpdat(lmpdat_b)
        sys_lmpdat_b.ts_coords.append(sys_dcd_ab.ts_coords[sys_lmpdat_a_natoms + 1:])
        sys_lmpdat_b.change_indices()
        sys_lmpdat_b.write_lmpdat(output_lmpdat_b, -1, title="Best frame of {} with index {}".format(os.path.basename(dcd_ab), index), cgcmm=True)
