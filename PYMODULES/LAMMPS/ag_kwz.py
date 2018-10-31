from __future__ import print_function, division
import os
#import copy
#import shutil as sl
import re
#import argparse
import math
#import time
import numpy as np
#from natsort import natsorted
#import itertools as it
#import scipy.stats
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
import pdb

#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator

#==============================================================================#
# Helper functions
#==============================================================================#


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

            finished_cycles = []
            folders_pwd = ["{}/{}".format(pwd, i) for i in os.listdir(pwd) if os.path.isdir(i)]

            # get last cycle from directory
            for folder in folders_pwd:
                cycle = re.match(r'.*?([0-9]+)$', folder).group(1)
                cycle = int(cycle)

                # avoid duplicates
                if cycle not in finished_cycles:
                    finished_cycles.append(cycle)

            return finished_cycles

        finished_cycles = get_finished_cycles()
        current_cycle = max(finished_cycles)
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
    total_cycles = range(total_cycles)
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


def berendsen_md(lmpcuts, group="all"):
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
    lmpcuts.berendsen(lmp, group=group)
    lmp.command("run {}".format(lmpcuts.runsteps))
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("reset_timestep 0")
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

    lmpcuts.nose_hoover(lmp, group=group)
    lmp.command("run {}".format(lmpcuts.runsteps))
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("reset_timestep 0")
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


def sysprep(lmpdat_out, lmpdat_main, lmpdat_add, dcd_add=None, frame_idx=-1):
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

    dcd_add : str (optional, default: None)
        dcd file to load on top of lmpdat_add

    frame_idx : int
        frame index of frame to add from dcd_add to lmpdat_add

    Returns
    -------
    success : bool
        True if successful, False otherwise.

    Writes a new lammps data file called sysprep_out_index

    """
    # read and transpose the main sys to the origin
    main_sys = aglmp.read_lmpdat(lmpdat_main)
    main_sys.transpose_by_cog(-1, [0, 0, 0], copy=False)
    _natoms = len(main_sys.atoms)

    # read and transpose the add sys to the origin
    add_sys = aglmp.read_lmpdat(lmpdat_add, dcd_add, frame_idx)
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


################################################################################
# Quenching
################################################################################

def quench(lmpcuts, lmpdat_main):
    """
    """
    def get_natms():
        """
        """
        with open(lmpdat_main, "r") as f_in:
            line = f_in.readline()
            while line != '':
                line = f_in.readline()
                if "atoms" in line:
                    return int(line.split()[0])

    natoms_main_sys = get_natms()

    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp, neigh=True)

    lmp.file(lmpcuts.settings_file)
    # change box type to not be periodic
    lmp.command("boundary f f f")
    lmpcuts.load_system(lmp)
    lmp.command("velocity all create {} 8455461 mom yes rot yes dist gaussian".format(lmpcuts.tstart))
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp)
    #lmp.command("variable n_hbond equal c_hb[1]")
    #lmp.command("variable E_hbond equal c_hb[2]")
    #lmp.command("thermo_style custom step temp epair v_E_hbond")

    # thermo output
    # compute dreiding hbonds and energies, if the pair style dreiding is set
    thermo_style = "thermo_style custom" + " ".join(lmpcuts.thermargs)

    thermo_dreiding = (
        "if $(is_active(pair_style,hbond/dreiding/lj)) then " +
        "\"" +
        "compute hb all pair hbond/dreiding/lj\n" +
        "variable n_hbond equal c_hb[1]\n" +
        "variable E_hbond equal c_hb[2]\n" +
        thermo_style +
        "v_n_hbond" + "v_E_hbond" +
        "\"" +
        "else " +
        "\"" +
        thermo_style +
        "\"")

    lmp.command(thermo_dreiding)
    pdb.set_trace()

    lmp.command("compute hb all pair hbond/dreiding/lj")
    lmpcuts.thermo(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

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
        prep_sys = aglmp.read_lmpdat(lmpcuts.input_lmpdat)
        cog = agm.get_cog(prep_sys.ts_coords[-1][natoms_main_sys + 1:])
        cog /= np.linalg.norm(cog, axis=0)  # unit vector
        cog *= -1
    else:
        cog = None

    cog = comm.bcast(cog, 0)
    # barostatting, thermostatting
    lmp.command("fix integrator grp_add_sys nvt temp {0} {1} 0.1".format(lmpcuts.tstart, lmpcuts.tstop))
    quench_success = False

    # 20 attempts to dock the molecule
    for _ in xrange(20):
        lmp.command(("fix force grp_add_sys addforce {0} {1} {2} every 100000").format(*cog))
        lmp.command("run {}".format(lmpcuts.runsteps))
        # remove pushing force
        lmp.command("unfix force")

        # post-optimization
        lmp.command("min_style quickmin")
        lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

        # check aggregate
        if rank == 0:
            quench_sys = aglmp.read_lmpdat(lmpcuts.input_lmpdat, dcd=lmpcuts.output_dcd)
            quench_success = quench_sys.check_aggregate()

        quench_success = comm.bcast(quench_success, 0)

        if quench_success is True:
            # write restart file
            lmpcuts.unfix_undump(pylmp, lmp)
            lmp.command("reset_timestep 0")
            lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
            lmp.command("clear")
            lmp.close()
            break

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


def _check_imaginary_clashes(sys_both, sys_a, sys_b, dcd_b):
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
    solvate_sys = aglmp.read_lmpdat(lmpdat_solvate, dcd_solvate)
    radii_mol, cogs_mol = _molecules_radii(solvate_sys)
    indent_strs = _fix_indent_ids(radii_mol, cogs_mol, "molecule", scale_start=10, scale_stop=15)

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
    _lmp_indent(lmp, indent_strs, lmpcuts.runsteps, keep_last_fixes=True)

    # solution system
    solvent_sys = aglmp.read_lmpdat(lmpcuts.input_lmpdat)
    solution_sys = mdu.merge_systems([solvate_sys, solvent_sys])
    solution_sys.reset_cells()

    # check solvate atoms with too close contacts to solvent atoms
    close_atoms = _check_imaginary_clashes(solution_sys, solvate_sys, solvent_sys, lmpcuts.output_dcd)
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
            indent_strs = _fix_indent_ids(all_radii_atoms, all_cogs_atoms, "atom", scale_start=factor_start, scale_stop=factor_stop)
            _lmp_indent(lmp, indent_strs, lmpcuts.runsteps, keep_last_fixes=False)
            close_atoms = _check_imaginary_clashes(solution_sys, solvate_sys, solvent_sys, lmpcuts.output_dcd)

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

    return close_atoms == []


def requench(lmpcuts):
    """
    """
    lmp = lammps()
    pylmp = PyLammps(ptr=lammps)
    lmp.file(lmpcuts.settings_file)
    # read restart file from previous run
    lmpcuts.read_system(lmp)
    lmpcuts.thermo(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    lmpcuts.minimize(lmp, style="cg", box_relax=False)
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()

def _append_data(data, lmplog):
    """
    """
    cur_log = agl.LmpLog()
    cur_log.read_lmplog(lmplog)
    cur_pes = cur_log.data[-1]["c_pe"]
    data.extend(cur_pes)


def vmd_rmsd_and_cluster(lmpdat_in, dcd_files, atm_idxs, xyz_out, percentage_to_check=80):
    """
    Align solvate molecules and calculate the rmsd.

    The first frame serves as the reference frame for the alignment and the rmsd
    calculation. The vmd modules molecule and atomsel have to be loaded for
    this function to work (and the coordinates had also to be read already).

    Parameters
    ----------
    lmpdat_in : str
        name of lammps data file with topology to load
    dcd_files : list of str
        file names of the dcd files to load into vmd

    Returns
    -------
    rmsds : list of floats
        rmsd values for each aligned selection
    cluster0 : list of ints
        indices of frames from the most populated cluster

    """
    import vmd
    import molecule
    import atomsel

    percentage_to_check /= 100.0
    # read topology and coordinates
    molecule.load("lmpdat", lmpdat_in)

    for dcd_file in dcd_files:
        molecule.read(0, "dcd", dcd_file, waitfor=-1)

    # rmsd with alignment for the solvent atoms
    rmsds = []
    selection = "index {}".format(" ".join(map(str, atm_idxs)))
    sel1 = atomsel.atomsel(selection, frame=0)

    frame_first = int(1 - percentage_to_check * molecule.numframes(0))
    frame_last = molecule.numframes(0)

    for frame_idx in range(frame_first, frame_last):
        sel2 = atomsel.atomsel(selection, frame=frame_idx)
        matrix_algin = sel2.fit(sel1)
        sel2.move(matrix_algin)
        rmsd = sel2.rmsd(sel1)
        rmsds.append(rmsd)

    # clustering (currently only possible using the tcl interface)
    tcl_selection = 'set tcl_sel1 [atomselect {} "index 0 to {} frame {}"]'
    tcl_sel1 = tcl_selection.format(0, num_solvate_atoms - 1, 0)
    vmd.evaltcl(tcl_sel1)
    cluster = "measure cluster $tcl_sel1 num 5 distfunc rmsd cutoff 1.0 first {} last {} step 1 selupdate True"
    clusters = vmd.evaltcl(cluster.format(frame_first, frame_last))
    clusters = re.findall("\{([^[\}]*)\}", clusters)
    cluster_0 = [int(i) for i in clusters[0].split()]
    molecule.write(0, "xyz", xyz_out, beg=cluster_0, end=cluster_0, selection=sel1)
    vmd.VMDexit("closing vmd")
    return (rmsds, cluster_0)


################################################################################
# Annealing
################################################################################

def _anneal(lmpcuts, pe_atm_idxs):
    """
    Use nose-hoover baro- and thermostat for productive annealing run.
    """

    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))

    # caveat: only possible with neigh no if calculating pe/atom on the graphics
    # card
    if lmpcuts.gpu is True:
        lmpcuts.use_gpu(lmp, neigh=False)

    lmp.file(lmpcuts.settings_file)
    lmpcuts.read_system(lmp)
    lmpcuts.thermargs.append("c_pe_solvate")
    lmpcuts.thermo(lmp)
    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
    lmpcuts.dump(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    # compute potential energy of the solvate (pair, bond, ...)
    # in order to compute the pe of the solvent, the solvent atom-ids must be
    # known
    pe_atm_ids = [i + 1 for i in range(pe_atm_idxs)]
    atm_ids_str = " ".join(map(str, pe_atm_ids))
    lmp.command("group resname_atoms " + atm_ids_str)
    lmp.command("compute pe_solvate resname_atoms pe/atom")
    lmp.command("compute pe resname_atoms reduce sum c_pe_solvate")
    lmpcuts.nose_hoover(lmp)
    lmp.command("run {}".format(lmpcuts.runsteps))
    lmpcuts.unfix_undump(pylmp, lmp)
    lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")
    lmp.close()

def _test_anneal_equil(data):
    """
    Check normality distribution of the underlying data.

    If the qq-plot fails, the distribution is not normal, if the other tests
    fail, normality can still be accepted if the qq plot states so.

    Parameters
    ----------
    data : list of floats

    Returns
    -------
    bool

    """
    qq_normal = ags.test_gauss_shape("qq", data)
    skew_normal = ags.test_gauss_shape("skewness", data)

    try:
        kurtosis_normal = ags.test_gauss_shape("kurtosis", data)
    except Warning:
        pass

    return qq_normal and (skew_normal or kurtosis_normal)

def anneal_productive(lmpcuts, atm_idxs_solvate, percentage_to_check):
    """
    #TODO check pressure, temperature equilibration?
    #TODO print picture of rmsd w/ cluster coloring?
    """
    #all_rmsds = []
    all_pe = []
    solvate_sys_natoms = len(atm_idxs_solvate)
    dcd_files = []
    log_files = []

    for run_idx in xrange(5):
        # rename current _anneal attempt
        lmpcuts.output_dcd = "{}_{}".format(run_idx, lmpcuts.output_dcd)
        lmpcuts.output_lmplog = "{}_{}".format(run_idx, lmpcuts.output_lmplog)

        # read previous log- and dcd-files if they exist already
        if os.path.isfile(lmpcuts.output_dcd) and os.path.isfile(lmpcuts.output_lmplog):
            dcd_files.append(lmpcuts.output_dcd)
            log_files.append(lmpcuts.output_lmplog)
            # read potential energy of the solvate and all coordinates
            if rank == 0:
                _append_data(all_pe, lmpcuts.output_lmplog)
                #molecule.read(0, "dcd", lmpcuts.output_dcd, waitfor=-1)
            continue

        else:
            _anneal(lmpcuts, atm_idxs_solvate)
            dcd_files.append(lmpcuts.output_dcd)
            log_files.append(lmpcuts.output_lmplog)

        if rank == 0:
            _append_data(all_pe, lmpcuts.output_lmplog)
            num_frames_to_check = percentage_to_check / 100 * len(all_pe)
            pe_normal = _test_anneal_equil(all_pe[num_frames_to_check:])

            # check if aggregate is still fine after the annealing run, i.e. last
            # frame has a solvate aggregate that is still complete
            aggregate_ok = False

            if pe_normal is True:
                solution_sys = aglmp.read_lmpdat(lmpcuts.input_lmpdat, lmpcuts.output_dcd)
                solution_sys_atoms_idxs = range(len(solution_sys.atoms))
                aggregate_ok = solution_sys.check_aggregate(excluded_atm_idxs=solution_sys_atoms_idxs[solvate_sys_natoms:])
        else:
            aggregate_ok = False

        aggregate_ok = comm.bcast(aggregate_ok)

        if aggregate_ok:
            break

    return (aggregate_ok, dcd_files, log_files)
