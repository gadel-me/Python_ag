#!/usr/bin/env python
from __future__ import print_function, division
from mpi4py import MPI
import os
import copy
import shutil as sl
import re
import argparse
import math
import numpy as np
from natsort import natsorted
#import itertools as it
import scipy.stats
from lammps import lammps, PyLammps
import Transformations as cgt
import md_elements as mde
import md_box as mdb
import ag_unify_md as agum
import ag_geometry as agm
import ag_unify_log as agul
import time

"""
Kawska Zahn Approach with lammps. Do not use neigh yes since it leads to segment-
ation faults. Always clear lammps or this will also lead to segmentation faults.
Clearing lammps not necessary if running solely on the cpu.
"""

"""
CAVEAT: DO NOT UNWRAP SOLVENT BOX, IT MUST BE WRAPPED OR OTHERWISE THE DENSITY
        IS WRONG DUE TO THE PERIODIC BOUNDARY CONDITIONS (MOLECULES THAT LEAVE
        THE BOX ARE OUTSIDE AND DO NOT GET MIRRORED BACK -> SOLVATE SHOULD BE IN
        THE MIDDLE OF THE BOX TO NOT BE WRAPPED AS WELL!)

CAVEAT: THE PATTERN STUFF IS ADDED SHOULD BE ACCORDING TO A CERTAIN PROBABILITY.
"""

#TODO When restarting, check if last run equilibrated successfully

#==============================================================================#
# Helper functions and variables
#==============================================================================#
_pwd = os.getcwd()
# percent of last values from log file to check
percentage_to_check = 80
_sqrt_3 = math.sqrt(3)

thermargs = ["step", "temp", "press", "vol", "density",
             "cella", "cellb", "cellc", "cellalpha", "cellbeta", "cellgamma",
             "etotal", "pe", "evdwl", "ecoul", "ebond", "eangle", "edihed", "eimp",
             "enthalpy"]


def create_folder(folder):
    """
    Create folder or skip creation if it already exists
    """
    try:
        os.mkdir(folder)
    except OSError:
        print("***Info: Folder {} already exists!".format(folder))


def write_to_log(s, filename="kwz_log"):
    """
    Write string 's' to log file named 'kwz.log'.
    """
    with open("kwz_log", "a") as kwz_log:
        kwz_log.write(s)


def process_data(data_lammps, dcd_lammps=None):
    """
    Read data files and merge, if solvent data and relaxed solvent data were given.

    Input:
        >   data_prepout            str; lammps data file (output from sys-preparation)
        >   data_solv               str; lammps data file (solvent)
        >   dcd_solvent_voids       str; dcd-file (output from creating voids into lammps data file)
    """
    # read output from quenching
    system = agum.Unification()
    system.read_lmpdat(data_lammps)

    if dcd_lammps is not None:
        # clean coordinates and boxes (saves memory)
        system.ts_coords = []
        system.ts_boxes  = []
        system.import_dcd(dcd_lammps)
        system.read_frames(frame=-2, to_frame=-1)  # read last frame
        system.close_dcd()

    return system


def check_energy_convergence(logfiles,
                             keyword="PotEng",
                             perecentage=100,
                             debug=False):
    """
    Calculate the skewness, p-value and z-score, to check if the values
    of 'keyword' are normally distributed (i.e. the MD-Simulation run ended
    successfully). "Generally speaking, the p-value is the probability of an
    outcome different that what was expected given the null hypothesis - in this
    case the probability of getting a skewness different from that of a normal
    distribution (which is 0 because of symmetry).

    Input:
        > logfiles      list; all written lammps-logfiles
        > keyword       str; keywords from thermo-output
        > perecentage   int; number of last values in %,
                        e.g. 80 means last 80 % of all values
        > debug         boolean; enable debug messaging
        > stage         str; which stage of the kwz-approach (relevant for
                        debugging)
    Output:
        > skewness, pvalue, zscore
        > data          list; all data for keyword from all files given
    """
    if debug:
        print("***check_energy_convergence", logfiles)

    normal_distribution = False
    data = []

    # gather all values from all logfiles given
    log_data = agul.LogUnification()
    log_data.read_lmplog(*logfiles)
    data = [value for data_index in xrange(len(log_data.data)) for value in
            log_data.data[data_index][keyword]]
    num_values = len(data)
    perecentage /= 100  # perecentage to per cent
    testdata = data[-int(perecentage*num_values):]
    #print(testdata)
    min_value_testdata = min(testdata)

    if debug is True:
        print(("***Info: Smallest value of "
               "{} % of all measured data is: {}").format(perecentage*100,
                                                          min_value_testdata)
              )

    skewness = scipy.stats.skew(testdata)
    zscore_skewness, pvalue_skewness = scipy.stats.skewtest(testdata)

    # normality test
    statistic_normality, pvalue_normality = scipy.stats.normaltest(testdata)

    if debug is True:
        print(("***Info: "
               "Skewness {}, "
               "Pvalue (skewness) {}, "
               "Pvalue (normality) {}").format(skewness,
                                               pvalue_skewness,
                                               pvalue_normality))

    # skewness ~ 0 with p-value <= 0.5 (normal distribution achieved)
    if (skewness < 0.3 and skewness > -0.3) and pvalue_skewness <= 0.10 and pvalue_normality <= 0.10:
        normal_distribution = True

    return (skewness, pvalue_skewness, pvalue_normality,
            normal_distribution, data, min_value_testdata)


def check_aggregate(mdsys, frame_id=-1, atm_atm_dist=4, unwrap=False, debug=False):
    """
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
    cube_side = atm_atm_dist/_sqrt_3  # -> atm_atm_dist/(2*math.sqrt(3)) * 2

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


#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator

#==============================================================================#
# Argument Parsing
#==============================================================================#

parser = argparse.ArgumentParser(
    prog="kawska_zahn.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Kawska-Zahn-Approach for accelerated crystallization simulations.")

# exclusive arguments (deprecated?)
#start = parser.add_mutually_exclusive_group()

#==================#
# general arguments
#==================#

parser.add_argument("-lmpm",
                    default=None,
                    metavar="*.lmpdat",
                    help="Lammps' data-file of the main system."
                    )

parser.add_argument("-lmpa",
                    default=None,
                    nargs="*",
                    metavar="*.lmpdat",
                    help="Lammps' data-file with atom-cube_sidetypes and coordinates for one single molecule to add to the current system. Atom types have to be defined already by the first data/restart file loaded!"
                    )

parser.add_argument("-lmps",
                    default=None,
                    metavar="*.lmpdat",
                    help="Create a solvent box for MD-Simulation. Data-file with one single solvent molecule to add."
                    )

parser.add_argument("-set",
                    metavar="*.lmpcfg",
                    required=True,
                    help="lammps' input-file/-script with basic simulation " +
                         "settings"
                    )

parser.add_argument("-pair_coeffs",
                    default=None,
                    metavar="*.lmpcfg",
                    required=True,
                    help="lammps'  script with lj, dreiding, etc. parameters"
                    )

parser.add_argument("-logsteps",
                    type=int,
                    default=1000,
                    help="log thermodynamic-, steps- and restart-files every" +
                         "'logsteps' steps"
                    )

parser.add_argument("-gpu",
                    default=False,
                    action="store_true",
                    help="utilize lammps' GPU package.",
                    )

parser.add_argument("-cycles",
                    type=int,
                    default=5,
                    help="Number of aggregations.")

parser.add_argument("-timeout",
                    metavar="00:01:00",
                    default="00:00:05",
                    help="allowed duration of simulation, for resubmitting purposes;  should be < 24h")

parser.add_argument("-pa",
                    "-pattern_add",
                    nargs="*",
                    default=0,
                    type=int,
                    help="The pattern in which looping order lmpa will be added, e.g. 0 1 2 3 3 1, repeats every 6 cycles")

#==========#
# quenching
#==========#

parser.add_argument("-quench_temp",
                    type=int,
                    default=5
                    )

parser.add_argument("-quench_steps",
                    type=int,
                    default=250000
                    )

parser.add_argument("-quench_logsteps",
                    type=int,
                    default=1000
                    )

#==========#
# annealing
#==========#

parser.add_argument("-void_solvent_steps",
                    type=int,
                    default=500
                    )

parser.add_argument("-start_relax_solvent_temp",
                    type=int,
                    default=10
                    )

parser.add_argument("-stop_relax_solvent_temp",
                    type=int,
                    default=10
                    )

parser.add_argument("-relax_solvent_steps",
                    type=int,
                    default=50000
                    )

parser.add_argument("-start_equil_anneal_temp",
                    type=int,
                    default=10
                    )

parser.add_argument("-stop_equil_anneal_temp",
                    type=int,
                    default=300,
                    help="Temperature at which system shall be annealed")

parser.add_argument("-equil_anneal_steps",
                    type=int,
                    default=500000
                    )

#parser.add_argument("-equil_anneal_ensemble",
#                    type=str,
#                    default="npt",
#                    help="nvt or npt")

parser.add_argument("-anneal_steps",
                    type=int,
                    default=2000000
                    )

parser.add_argument("-additional_anneal_steps",
                    type=int,
                    default=100000)

parser.add_argument("-anneal_logsteps",
                    type=int,
                    default=500)

parser.add_argument("-requench_steps",
                    type=int,
                    default=150000
                    )

args = parser.parse_args()

#==============================================================================#
# Remaining cycles and molecule to add pattern
#==============================================================================#

if rank == 0:
    finished_cycles  = []
    folders_pwd      = ["{}/{}".format(_pwd, i) for i in
                        os.listdir(_pwd) if os.path.isdir(i)]

    # get last cycle from directory
    for folder in folders_pwd:
        cycle = re.match(r'.*?([0-9]+)$', folder).group(1)
        cycle = int(cycle)
        # avoid duplicates
        if cycle not in finished_cycles:
            finished_cycles.append(cycle)

    del (folders_pwd)

    # requench_out does not exist -> current cycle has not finished yet
    if len(finished_cycles) >= 1:
        current_cycle = max(finished_cycles)
        requench_out = _pwd + "/requench_{}/".format(current_cycle) + "requench_{}.dcd".format(current_cycle)

        if os.path.isfile(requench_out) is True:
            next_cycle = current_cycle + 1
        else:
            next_cycle = current_cycle

        #del requench_out
    else:
        next_cycle = 0

    del (finished_cycles)

    total_cycles = range(args.cycles)
    idx_next_cycle = total_cycles.index(next_cycle)
    del (total_cycles)

    #===========================#
    # Molecule to add by pattern
    #===========================#

    id_pattern = 0
    num_patterns = len(args.pa) - 1  # indices start with 0 so one less
    total_cycles = []

    for cycle in range(args.cycles):
        if id_pattern > num_patterns:
            id_pattern = 0
        total_cycles.append((cycle, args.pa[id_pattern]))
        id_pattern += 1

    remaining_cycles = total_cycles[idx_next_cycle:]
    del (id_pattern, num_patterns, total_cycles, cycle)
else:
    remaining_cycles = None

remaining_cycles = comm.bcast(remaining_cycles, 0)

#==============================================================================#
# Kawska Zahn Approach
#==============================================================================#

for curcycle, idx_lmpa in remaining_cycles:

    #==========================================================#
    # Define folders and files, retrieve stage of current cycle
    #==========================================================#

    write_to_log("Cycle: {:d}\n".format(curcycle))

    # declare folder names for each cycle
    sysprep_dir  = _pwd + "/sysprep_{}/".format(curcycle)
    quench_dir   = _pwd + "/quench_{}/".format(curcycle)
    anneal_dir   = _pwd + "/anneal_{}/".format(curcycle)
    requench_dir = _pwd + "/requench_{}/".format(curcycle)

    # system preparation
    sysprep_out  = sysprep_dir + "sysprep_out_{}.lmpdat".format(curcycle)

    # quench
    quench_out = quench_dir + "quench_out_{}.lmprst".format(curcycle)
    quench_rst = quench_dir + "quench_rst_{}.lmprst".format(curcycle)
    quench_dcd = quench_dir + "quench_{}.dcd".format(curcycle)
    quench_log = quench_dir + "quench_{}.lmplog".format(curcycle)

    # anneal -> solvent
    void_solv_out = anneal_dir + "void_solv_{}".format(curcycle) + "_out.lmpdat"
    void_solv_rst = anneal_dir + "void_solv_{}".format(curcycle) + "_tmp.rst"
    void_solv_dcd = anneal_dir + "void_solv_{}".format(curcycle) + ".dcd"
    void_solv_log = anneal_dir + "void_solv_{}".format(curcycle) + ".lmplog"

    relax_solv_in  = anneal_dir + "relax_solv_{}".format(curcycle) + "_in.lmpdat"
    relax_solv_out = anneal_dir + "relax_solv_{}".format(curcycle) + "_out.lmprst"
    relax_solv_rst = anneal_dir + "relax_solv_{}".format(curcycle) + "_tmp.lmprst"
    relax_solv_dcd = anneal_dir + "relax_solv_{}".format(curcycle) + ".dcd"
    relax_solv_log = anneal_dir + "relax_solv_{}".format(curcycle) + ".lmplog"

    # anneal -> equilibration/heating
    equil_anneal_out = anneal_dir + "equil_anneal_{}".format(curcycle) + "_out.lmprst"
    equil_anneal_rst = anneal_dir + "equil_anneal_{}".format(curcycle) + "_tmp.lmprst"
    equil_anneal_dcd = anneal_dir + "equil_anneal_{}".format(curcycle) + ".dcd"
    equil_anneal_log = anneal_dir + "equil_anneal_{}".format(curcycle) + ".lmplog"

    # anneal -> productive
    solvate_anneal_out = anneal_dir + "anneal_{}".format(curcycle) + "_solvate_out.xyz"
    solvent_anneal_out = anneal_dir + "anneal_{}".format(curcycle) + "_solvent_out.xyz"
    anneal_rst = anneal_dir + "anneal_{}".format(curcycle) + "_tmp.lmprst"
    anneal_dat = anneal_dir + "anneal_{}".format(curcycle) + ".lmpdat"
    #anneal_dcd = anneal_dir + "anneal_{}".format(curcycle) + ".dcd"
    anneal_log = anneal_dir + "anneal_{}".format(curcycle) + ".lmplog"

    # requench
    tmp_solvate_anneal_out = requench_dir + "requench_{}".format(curcycle) + "_tmp_solvate_out.lmpdat"
    #requench_out = requench_dir + "requench_{}".format(curcycle) + "_out.lmpdat"
    requench_dcd = requench_dir + "requench_{}".format(curcycle) + ".dcd"
    requench_log = requench_dir + "requench_{}".format(curcycle) + ".lmplog"

    # important files from previous run
    pre_sysprep_out        = "{0}/sysprep_{1}/sysprep_out_{1}.lmpdat".format(_pwd, curcycle-1)
    pre_solvent_anneal_out = "{0}/anneal_{1}/anneal_{0}_solvent_out.xyz".format(_pwd, curcycle-1)
    pre_requench_dcd       = "{0}/requench_{1}/requench_{1}.dcd".format(_pwd, curcycle-1)

    if os.path.isfile(quench_out) is True:
        quench_success = True
    else:
        quench_success = False

    if os.path.isfile(solvate_anneal_out) is True:
        anneal_success = True
    else:
        anneal_success = False

    #==========================================================================#
    # Aggregation
    #==========================================================================#

    # define main system
    if os.path.isfile(pre_requench_dcd) is True:
        main_prep_sys = pre_sysprep_out
    else:
        main_prep_sys = os.path.abspath(args.lmpm)

    while anneal_success is False:
        while quench_success is False:

            #==================================================================#
            # 1. System Preparation
            #==================================================================#

            # very first run or last stage has finished successfully
            if os.path.isfile(sysprep_out) is False:

                if rank == 0:
                    create_folder(sysprep_dir)
                    write_to_log("\tStage: SYS-PREPARATION\n")

                    #TODO ADD NEW SYSTEM ACCORDING TO A CERTAIN PROBABILITY, NOT BY PATTERN
                    # define add system(s)
                    add_prep_sys = os.path.abspath(args.lmpa[idx_lmpa])

                    write_to_log("\tInput-File Main: {}\n".format(main_prep_sys))
                    write_to_log("\tInput-File Add: {}\n".format(add_prep_sys))

                    # add new molecule to the main system
                    for _ in xrange(100):
                        # process main system
                        main_sys = process_data(main_prep_sys)

                        if os.path.isfile(pre_requench_dcd) is True:
                            main_sys.ts_coords = []
                            main_sys.import_dcd(pre_requench_dcd)
                            main_sys.read_frames(frame=-2, to_frame=-1)
                            main_sys.close_dcd()

                        natoms_main_sys = len(main_sys.atoms)
                        main_sys.transpose_by_cog(0, [0, 0, 0], copy=False)  # relocate main system to origin
                        # process add system
                        add_sys = process_data(add_prep_sys)
                        add_sys.transpose_by_cog(0, [0, 0, 0], copy=False)  # relocate add system to origin (for proper placement on sphere)

                        # rotate add-system by an arbitrary rotation matrix,
                        # consisting of 3 arbitrary rotations A, B and C
                        M_arb_rot_x = agm.arb_rot_matrix([1, 0, 0])  # A
                        M_arb_rot_y = agm.arb_rot_matrix([0, 1, 0])  # B
                        M_arb_rot_z = agm.arb_rot_matrix([0, 0, 1])  # C
                        # multiply all three matrices; ABC=A(BC)=(AB)C
                        M_arb_xy  = np.matmul(M_arb_rot_x, M_arb_rot_y)  # AB
                        M_arb_xyz = np.matmul(M_arb_xy, M_arb_rot_z)  # ABC
                        # rotate the system
                        sys_add_natoms = len(add_sys.atoms)
                        sys_add_all_atms = range(sys_add_natoms)
                        add_sys.mm_atm_coords(0, M_arb_xyz, False, *sys_add_all_atms)
                        del (M_arb_rot_x, M_arb_rot_y, M_arb_rot_z, M_arb_xy, M_arb_xyz)

                        # transpose add system to an arbitrary point on the sphere
                        # radius of a sphere that enwraps the main system
                        main_sys_radius = main_sys.get_system_radius(0)

                        # radius of a sphere that enwraps the to-be-added system
                        add_sys_radius = add_sys.get_system_radius(0)

                        # get kawska-zahn-radius (kwz) (sum of both system's radii)
                        kwz_radius = main_sys_radius + add_sys_radius
                        kwz_radius += 1  # enlarge radius so atoms do not clash
                        rn_pos = agm.points_on_sphere(npoints=1, ndim=3, radius=kwz_radius)[0]
                        M_trans = cgt.translation_matrix(rn_pos)
                        add_sys.mm_atm_coords(0, M_trans, False, *sys_add_all_atms)
                        del (rn_pos, M_trans, sys_add_all_atms)

                        # merge both systems
                        main_sys.extend_universe(add_sys, mode="merge")
                        del add_sys

                        # add new orthogonal box
                        box_diameter  = 2 * (main_sys_radius*2 + add_sys_radius*2)
                        a = b = c = box_diameter

                        # DEBUGGING NOW
                        a = b = c = 100

                        alpha = beta = gamma = math.pi/2
                        main_sys.ts_boxes = []
                        main_sys.ts_boxes.append(mdb.Box(boxtype="lattice",
                                                 ltc_a=a, ltc_b=b, ltc_c=c,
                                                 ltc_alpha=alpha, ltc_beta=beta,
                                                 ltc_gamma=gamma))
                        del (main_sys_radius, add_sys_radius, box_diameter,
                             a, b, c, alpha, beta, gamma)
                        main_sys.ts_boxes[0].box_lat2lmp(triclinic=False)

                        # group atoms by bonds
                        main_sys.fetch_molecules_by_bonds()
                        main_sys.mols_to_grps()

                        # check interatomic distances
                        main_sys.create_linked_cells(-1)
                        close_atms = main_sys.chk_atm_dist(frame_id=-1, min_dist=2)

                        # replace new atoms if they clash with atoms from the main system
                        # (which should not happen, just in case we grow an ellipsoid)
                        close_atms = sorted(close_atms)

                        # check if new atoms were placed too near to the main system
                        try:
                            if close_atms[-1] > natoms_main_sys:
                                print("***Sysprep-Warning: New atoms were placed too near! Rebuilding...")
                            else:
                                break

                        except IndexError:  # close_atoms is empty array (-> no close atoms)
                            break

                        del (close_atms)
                    else:
                        print("***Sys-Prep-Error: Docking failed!")

                    del (natoms_main_sys)

                    # write new data file
                    main_sys.change_indices(incr=1, mode="increase")
                    main_sys.write_lmpdat(sysprep_out, frame_id=0, title="System ready" +
                                          "for docking", cgcmm=True)

                    # write successful steps to log file
                    write_to_log("\tSphere-radius: {}\n".format(kwz_radius))
                    write_to_log("\tOutput-File: {}\n\n".format(sysprep_out))

                    del (main_sys, kwz_radius)

            #===================================================================
            # 2. System Quenching
            #===================================================================

            if os.path.isfile(quench_out) is False:
                natoms_main_sys = None

                if rank == 0:
                    create_folder(quench_dir)
                    write_to_log("\tStage: QUENCHING\n")
                    # integrate only added atoms from sys-add
                    main_sys = process_data(main_prep_sys)
                    natoms_main_sys = len(main_sys.atoms)
                    del main_sys
                    out_quench_sys = process_data(sysprep_out)

                natoms_main_sys = comm.bcast(natoms_main_sys, 0)

                quench_lmp   = lammps()
                quench_pylmp = PyLammps(ptr=quench_lmp)
                quench_lmp.command("log {} append".format(quench_log))

                # load gpu package
                if args.gpu:
                    # neighbor list building on gpu (currently) leads to segmentation faults
                    #quench_lmp.command("package gpu 1 neigh yes")
                    quench_lmp.command("package gpu 1 neigh no")  # DEBUGGING
                    quench_lmp.command("suffix gpu")

                # load general settings (units, boundary, dimension, atom_style, etc. pp.)
                quench_lmp.file(args.set)

                # box needs not to be periodic here
                quench_lmp.command("boundary f f f")

                # define simulation temps
                quench_temp = args.quench_temp

                # read data or restart file
                if os.path.isfile(quench_rst) is True:
                    quench_lmp.command("read_restart {}".format(quench_rst))
                else:
                    quench_lmp.command("read_data {}".format(sysprep_out))
                    quench_lmp.command("velocity all create {} 483806 rot yes dist gaussian".format(quench_temp))

                if args.pair_coeffs is not None:
                    quench_lmp.file(args.pair_coeffs)

                quench_lmp.command("group loose id > {}".format(natoms_main_sys))

                # ice cube prevention
                quench_lmp.command("fix ic_prevention all momentum {} linear 1 1 1 angular rescale".format(100))

                # barostatting, thermostatting
                quench_lmp.command("fix ensemble_quench loose nvt temp {0} {0} 0.1".format(quench_temp))

                # set an additional push to the added atoms that should be docked
                # (prevents losing atoms due to being localized outside the cutoff)
                if rank == 0:
                    cog = agm.get_cog(out_quench_sys.ts_coords[-1][natoms_main_sys+1:])
                    cog   /= np.linalg.norm(cog, axis=0)  # unit vector
                    cog   *= -1e-5  # force vector of length 2 pointing towards the origin
                else:
                    cog = None

                del (natoms_main_sys)
                cog = comm.bcast(cog, 0)
                quench_lmp.command(("fix group_loose_addforce "
                                    "loose addforce {0} {1} {2}").format(*cog))

                del (quench_temp)

                # logging of thermodynamic data
                quench_lmp.command("thermo_style custom " + " ".join(thermargs))
                quench_lmp.command("thermo_modify lost warn flush yes")
                quench_lmp.command("thermo {}".format(args.quench_logsteps))

                # trajectory
                quench_lmp.command("dump QUENCH_DUMP all dcd {} {}".format(
                                   args.quench_logsteps*2, quench_dcd))
                # unwrap trajectory coordinates
                quench_lmp.command("dump_modify QUENCH_DUMP unwrap yes")

                # intermediate restart files
                #tmp1_rst  = quench_dir + "tmp-1-quench_out_{}.lmprst".format(curcycle)
                steps_rst = args.quench_logsteps * 10
                quench_lmp.command("restart {} {} {}".format(steps_rst, quench_rst, quench_rst))

                # initial quench run
                quench_steps = args.quench_steps
                # 1st run with a little push towards the aggregate
                if os.path.isfile(quench_rst) is False:
                    quench_lmp.command("run {}".format(int(quench_steps*0.2)))
                quench_lmp.command("unfix group_loose_addforce")
                # 2nd run without push bias
                quench_lmp.command("run {} upto".format(quench_steps))

                # check if additional molecule docked successfully
                if rank == 0:
                    out_quench_sys.import_dcd(quench_dcd)
                    out_quench_sys.read_frames(frame=-2, to_frame=-1)
                    out_quench_sys.close_dcd()
                    quench_success = check_aggregate(out_quench_sys,
                                                     frame_id=-1,
                                                     atm_atm_dist=4)
                    del out_quench_sys

                quench_success = comm.bcast(quench_success, root=0)

                ##=====================================#
                # quench energy check, additional runs
                ##=====================================#
                #total_additional_quench_runs = 5
                #chck_vals_pe_quench = int(quench_steps/args.quench_logsteps) - int(quench_steps/args.quench_logsteps*0.1)
#
                #for additional_quench_run in xrange(total_additional_quench_runs):
#
                #    # check if additional molecule docked successfully
                #    if rank == 0:
                #        out_quench_sys = process_data(sysprep_out, quench_dcd)
                #        quench_success = check_aggregate(out_quench_sys,
                #                                         atm_atm_dist=4,
                #                                         frame_id=0)
                #        del out_quench_sys
#
                #    quench_success = comm.bcast(quench_success, root=0)
#
                #    if quench_success is False:
                #        break
#
                #    # ENERGY CHECK NOT NECESSARY IF MOLECULE HAS DOCKED SUCCESSFULLY
                #    ## check quenching success by energy
                #    #normal_distributed_pe_quench = None
##
                #    #if rank == 0:
                #    #    (skewness_pe_quench,
                #    #     pvalue_skewness_pe_quench,
                #    #     pvalue_normality_quench,
                #    #     normal_distributed_pe_quench,
                #    #     _) = check_energy_convergence(quench_log,
                #    #                                   keyword="PotEng",
                #    #                                   stage="Quenching",
                #    #                                   percentage=chck_vals_pe_quench,
                #    #                                   debug=True)
                #    #    write_to_log(("Skewness {}, "
                #    #                  "Pvalue (skewness) {}, "
                #    #                  "Pvalue (normality) {}\n").format(skewness_pe_quench,
                #    #                                                    pvalue_skewness_pe_quench,
                #    #                                                    pvalue_normality_quench))
                #    #    del (skewness_pe_quench, pvalue_skewness_pe_quench, pvalue_normality_quench)
##
                #    #normal_distributed_pe_quench = comm.bcast(normal_distributed_pe_quench, 0)
##
                #    ## stop further runs, if normal distribution is achieved
                #    #if normal_distributed_pe_quench is True:
                #    #    break
                #    #else:
                #    #    if additional_quench_run < total_additional_quench_runs - 1:
                #    #        quench_lmp.command("run {}".format(quench_steps))
#
                ##del (quench_steps, chck_vals_pe_quench)

                # write data/restart file, adjust lammps written data file,
                # rename lammps log file
                if quench_success is True:
                    if rank == 0:
                        print("***Quench-Info: Molecules successfully docked.")

                    # write restart file
                    quench_lmp.command("undump QUENCH_DUMP")
                    quench_lmp.command("unfix ensemble_quench")
                    quench_lmp.command("unfix ic_prevention")
                    quench_lmp.command("reset_timestep 0")
                    quench_lmp.command("write_restart {}".format(quench_out))
                else:
                    sysprep_success = False

                    if rank == 0:
                        print("***Quench-Warning: System did not aggregate! Quenching failed! Removing folders!")
                        sl.rmtree(sysprep_dir)
                        sl.rmtree(quench_dir)
                        write_to_log("\tDocking failed\n")

                # delete unnecessary files
                if rank == 0:
                    try:
                        os.remove(quench_rst)
                    except OSError:
                        pass

                quench_lmp.command("clear")
                quench_lmp.close()
                del quench_lmp

                if rank == 0:
                    print("***Quenching-Info: Quenching done!")

        #======================================================================#
        # 3. ANNEALING
        #======================================================================#

        if os.path.isfile(solvate_anneal_out) is False:
            if rank == 0:
                create_folder(anneal_dir)
                write_to_log("\tStage: ANNEALING\n")

            # solvate
            if rank == 0:
                solvate_sys = process_data(sysprep_out, quench_dcd)
                natoms_solvate_sys = len(solvate_sys.atoms)
            else:
                natoms_solvate_sys = None

            natoms_solvate_sys = comm.bcast(natoms_solvate_sys, 0)

            if args.lmps is None:
                if rank == 0:
                    # solution is solvate + vacuum
                    solution_sys = solvate_sys
            else:

                #================#
                # Prepare solvent
                #================#

                # solvent
                if rank == 0:
                    solvent_sys = process_data(args.lmps)
                    natoms_solvent_sys = len(solvent_sys.atoms)

                    # read geometry from previous annealing
                    if os.path.isfile(pre_solvent_anneal_out) is True:
                        solvent_sys.ts_coords = []
                        solvent_sys.read_xyz(pre_solvent_anneal_out)

                    # solute + solvent
                    solution_sys = copy.deepcopy(solvate_sys)
                    solution_sys.extend_universe(solvent_sys,
                                                 u1_frame_id=-1,
                                                 u2_frame_id=-1,
                                                 mode="merge")
                    solution_sys.ts_boxes  = []
                    solution_sys.ts_coords = []
                    solution_sys.mix_pair_types(mode="ii",
                                                mix_style="arithmetic")
                    solution_sys.fetch_molecules_by_bonds()
                    solution_sys.mols_to_grps()

                #========================#
                # solvent void generation
                #========================#

                if os.path.isfile(void_solv_out) is False:

                    # get center of geometry and sphere of each solvate molecule
                    radii_sphere = []
                    cogs_solvate = []

                    if rank == 0:
                        for molecule in solvate_sys.molecules:
                            sphere, cog = solvate_sys.get_mol_radius(-1, *molecule)
                            radii_sphere.append(sphere)
                            cogs_solvate.append(cog)

                        del (sphere, cog, molecule)
                        num_solvate_molecules  = len(solvate_sys.molecules)
                    else:
                        num_solvate_molecules = None

                    num_solvate_molecules = comm.bcast(num_solvate_molecules, 0)
                    radii_sphere = comm.bcast(radii_sphere, 0)
                    cogs_solvate = comm.bcast(cogs_solvate, 0)

                    # start lammps
                    void_lmp   = lammps()
                    void_pylmp = PyLammps(ptr=void_lmp)
                    void_lmp.command("log {} append".format(void_solv_log))

                    if args.gpu is True:
                        #void_lmp.command("package gpu 1 neigh yes")
                        void_lmp.command("package gpu 1 neigh no")  # DEBUGGING
                        void_lmp.command("suffix gpu")

                    # read settings file
                    void_lmp.file(args.set)

                    # load last void state
                    if os.path.isfile(void_solv_rst) is True:  # void run was aborted, restart file written; start from here
                        void_lmp.command("read_restart {}".format(void_solv_rst))
                    elif os.path.isfile(pre_solvent_anneal_out) is True:  # output from previous run
                        void_lmp.command("read_data {}".format(void_solv_out))
                    else:  # first run with solvent
                        void_lmp.command("read_data {}".format(args.lmps))

                    if args.pair_coeffs is not None:
                        void_lmp.file(args.pair_coeffs)

                    # load further lammps settings
                    void_lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")
                    void_lmp.command("thermo_style custom " + " ".join(thermargs))
                    void_lmp.command("thermo_modify lost warn flush yes")
                    void_lmp.command("thermo {}".format(args.logsteps))

                    # write restart, dcd files
                    void_lmp.command("restart {0} {1} {1}".format(args.logsteps*50, void_solv_rst))
                    void_lmp.command("dump void_dcd all dcd {} {}".format(args.logsteps, void_solv_dcd))
                    void_lmp.command("dump_modify void_dcd unwrap yes")
                    void_lmp.command(("fix create_void {} npt temp {} {} 0.1 "
                                      "iso 1.0 1.0 1").format("all", 10, 10))

                    void_steps = args.void_solvent_steps
                    # create magnifying voids (reduces risk of entrapped solvent molecules)
                    molecules_names_void_fix = ["void_indent_molecule_{}".format(index_void_fix) for
                                                index_void_fix in xrange(num_solvate_molecules)]
                    del (index_void_fix, num_solvate_molecules)

                    #=============================================#
                    # initial spheres around all solvate molecules
                    #=============================================#

                    for molecule_sphere_run in xrange(1, 12):

                        # set fixes for void_indent
                        for molecule_void_fix_name, molecule_void_sphere, molecule_void_cog in zip(molecules_names_void_fix,
                                                                                                   radii_sphere,
                                                                                                   cogs_solvate):
                            growing_molecule_sphere = (molecule_void_sphere *
                                                       0.1 * molecule_sphere_run)
                            molecule_void_fix = "fix {0} all indent 1 sphere {c[0]} {c[1]} {c[2]} {1} side out".format(
                                molecule_void_fix_name,
                                growing_molecule_sphere,
                                c=molecule_void_cog)
                            void_lmp.command(molecule_void_fix)

                        void_lmp.command("run {}".format(void_steps))

                        # delete all fixes in order to alter them in the next run
                        # do not delete them after the last run
                        if molecule_sphere_run < 11:
                            for molecule_void_fix_name in molecules_names_void_fix:
                                void_lmp.command("unfix {0}".format(molecule_void_fix_name))

                        del (molecule_void_fix_name, molecule_void_sphere,
                             molecule_void_cog, growing_molecule_sphere)

                    del molecule_sphere_run

                    #==============================================#
                    # additional spheres around close solvate atoms
                    #==============================================#

                    close_atoms_solvate = []
                    addnames_void_fix = []

                    atoms_coords_void = []
                    atoms_radii_void  = []
                    atoms_fixes_void  = []

                    # 50 attempts to make stuff grow large enough
                    for _ in xrange(50):

                        #================================================#
                        # interatomic distances (solvate - solvent) check
                        #================================================#

                        # array to test if any atoms still have close contacts
                        current_close_atoms_solvate = []

                        if rank == 0:
                            solution_sys.ts_coords  = []
                            solution_sys.ts_boxes   = []
                            solution_sys.ts_lnk_cls = []

                            # read latest solvent coordinates and boxes
                            solvent_sys.import_dcd(void_solv_dcd)
                            solvent_sys.read_frames(frame=-2, to_frame=-1)
                            solvent_sys.close_dcd()

                            # concatenate solute and latest solvent coordinates
                            solution_sys.ts_coords.append(
                                np.concatenate((solvate_sys.ts_coords[-1], solvent_sys.ts_coords[-1]))
                            )
                            solution_sys.ts_boxes = solvent_sys.ts_boxes
                            solution_sys.create_linked_cells(-1,
                                                             rcut_a=2,
                                                             rcut_b=2,
                                                             rcut_c=2)
                            close_atoms_solution = solution_sys.chk_atm_dist(-1,
                                                                             min_dist=2.0,
                                                                             exclude_same_molecule=True)
                            solution_sys.ts_lnk_cls = []
                            solution_sys.ts_coords  = []
                            solution_sys.ts_boxes   = []

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

                                print(void_pylmp.fixes)
                                time.sleep(5)

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

                    # write final coordinates and topology as lammps data file
                    if rank == 0:
                        solvent_sys.import_dcd(void_solv_dcd)
                        solvent_sys.read_frames(frame=-2, to_frame=-1)
                        solvent_sys.close_dcd()
                        # append solvate and updated solvent coordinates to the solution
                        solution_sys.ts_coords.append(
                            np.concatenate((solvate_sys.ts_coords[-1],
                                            solvent_sys.ts_coords[-1]))
                        )
                        # adapt box from solvent
                        solution_sys.ts_boxes = solvent_sys.ts_boxes
                        solution_sys.change_indices(incr=1, mode="increase")
                        solution_sys.ts_boxes[-1].lmp_xy = None
                        solution_sys.ts_boxes[-1].lmp_xz = None
                        solution_sys.ts_boxes[-1].lmp_yz = None
                        solution_sys.write_lmpdat(void_solv_out, -1, cgcmm=True,
                                                  title="Solvent-Voids with solvate, ready for solvent relaxation.")
                        solution_sys.change_indices(incr=1, mode="decrease")

                #===================#
                # solvent relaxation
                #===================#

                if os.path.isfile(relax_solv_out) is False:

                    if rank == 0:
                        write_to_log("\tSolvent-Relaxation\n")

                    relax_solvent_lmp = lammps()
                    relax_solvent_lmp.command("log {} append".format(relax_solv_log))

                    # load gpu package
                    if args.gpu:
                        #relax_solvent_lmp.command("package gpu 1 neigh yes")
                        relax_solvent_lmp.command("package gpu 1 neigh no")
                        relax_solvent_lmp.command("suffix gpu")

                    relax_solvent_lmp.file(args.set)

                    start_relax_temp = args.start_relax_solvent_temp
                    stop_relax_temp  = args.stop_relax_solvent_temp

                    if os.path.isfile(relax_solv_rst) is True:
                        #sl.copy2(relax_solv_dcd, relax_solv_dcd + "_pre_run.dcd")
                        relax_solvent_lmp.command("read_restart {}".format(relax_solv_rst))
                    else:
                        relax_solvent_lmp.command("read_data {}".format(void_solv_out))
                        relax_solvent_lmp.command(("velocity all create {} 483806 "
                                                   "rot yes dist gaussian").format(start_relax_temp))

                    if args.pair_coeffs is not None:
                        relax_solvent_lmp.file(args.pair_coeffs)

                    # read last temperature  # DEBUGGING not necessary since
                    # this is to relax the solvent before annealing begins
                    #if os.path.isfile(relax_solv_rst) is True:
                    #    relax_log = agul.LogUnification()
                    #    relax_log.read_lmplog(relax_log)
                    #    # reassign starting temperature due to last run
                    #    start_relax_temp = relax_log.data[-1]["Temp"][-1]

                    relax_solvent_lmp.command("thermo_style custom " + " ".join(thermargs))
                    relax_solvent_lmp.command("thermo_modify lost warn flush yes")
                    relax_solvent_lmp.command("thermo {}".format(args.logsteps))
                    relax_solvent_lmp.command("dump dump_relax_solvent all dcd {} {}".format(args.logsteps, relax_solv_dcd))
                    relax_solvent_lmp.command("dump_modify dump_relax_solvent unwrap yes")
                    relax_solvent_lmp.command(("fix ic_prevention all momentum {} "
                                               "linear 1 1 1 angular rescale").format(100))
                    # lammps indices start with 1 -> '<=' number atoms of solvate sys
                    relax_solvent_lmp.command("group solvate id <= {}".format(natoms_solvate_sys))
                    relax_solvent_lmp.command("group solvent id > {}".format(natoms_solvate_sys))
                    relax_solvent_lmp.command("restart {} {} {}".format(args.logsteps*25,
                                                                        relax_solv_rst,
                                                                        relax_solv_rst))

                    # berendsen thermo- and barostatting
                    #relax_solvent_lmp.command(("fix fix_nve_relax_solvent {} nve").format("solvent"))
                    #relax_solvent_lmp.command(("fix fix_temp_berendsen_relax_solvent {} "
                    #                           "temp/berendsen {} {} 0.1").format("solvent",
                    #                          start_relax_temp, stop_relax_temp))
                    #relax_solvent_lmp.command(("fix fix_press_berendsen_relax_solvent {} "
                    #                           "press/berendsen iso 1.0 1.0 1.0").format("solvent"))

                    # hoover-nose thermo- and barostatting
                    relax_solvent_lmp.command(("fix ensemble_relax_solvent {} "
                                               "npt temp {} {} 0.1 "
                                               "iso 1.0 1.0 1").format("solvent",
                                                                       start_relax_temp,
                                                                       stop_relax_temp))

                    relaxation_steps = args.relax_solvent_steps
                    relax_solvent_lmp.command("run {} upto".format(relaxation_steps))
                    del relaxation_steps
                    relax_solvent_lmp.command("undump dump_relax_solvent")
                    relax_solvent_lmp.command("unfix ensemble_relax_solvent")
                    #relax_solvent_lmp.command("unfix fix_nve_relax_solvent")
                    #relax_solvent_lmp.command("unfix fix_temp_berendsen_relax_solvent")
                    #relax_solvent_lmp.command("unfix fix_press_berendsen_relax_solvent")
                    relax_solvent_lmp.command("unfix ic_prevention")
                    relax_solvent_lmp.command("reset_timestep 0")
                    relax_solvent_lmp.command("write_restart {}".format(relax_solv_out))
                    relax_solvent_lmp.command("clear")
                    relax_solvent_lmp.close()

            #=========================================#
            # annealing (solution) - equilibration run
            #=========================================#

            start_anneal_temp = args.start_equil_anneal_temp
            stop_anneal_temp  = args.stop_equil_anneal_temp

            # read last temperature
            if os.path.isfile(equil_anneal_rst) is True and os.path.isfile(equil_anneal_log) is True:
                eq_log = agul.LogUnification()
                eq_log.read_lmplog(equil_anneal_log)
                # reassign starting temperature due to last run
                start_anneal_temp = eq_log.data[-1]["Temp"][-1]

            equil_steps = args.equil_anneal_steps

            if os.path.isfile(equil_anneal_out) is False:
                equil_anneal_lmp = lammps()
                equil_anneal_lmp.command("log {} append".format(equil_anneal_log))

                # load gpu package
                if args.gpu:
                    #equil_anneal_lmp.command("package gpu 1 neigh yes")
                    equil_anneal_lmp.command("package gpu 1 neigh no")
                    equil_anneal_lmp.command("suffix gpu")

                # load settings (units, boundary, dimension, atom_style, etc. pp.)
                equil_anneal_lmp.file(args.set)

                # in vacuo -> no boundary conditions
                if args.lmps is None:
                    equil_anneal_lmp.command("boundary f f f")

                # in solvent or in vacuo; restart
                if os.path.isfile(equil_anneal_rst) is True:
                    equil_anneal_lmp.command("read_restart {}".format(equil_anneal_rst))
                # in solvent; 1st run
                elif args.lmps is not None:
                    equil_anneal_lmp.command("read_restart {}".format(relax_solv_out))
                    equil_anneal_lmp.command("velocity all create {} 483806 rot yes dist gaussian".format(start_anneal_temp))
                # in vacuo; 1st run
                else:
                    equil_anneal_lmp.command("read_restart {}".format(quench_out))
                    equil_anneal_lmp.command("velocity all create {} 483806 rot yes dist gaussian".format(start_anneal_temp))

                if args.pair_coeffs is not None:
                    equil_anneal_lmp.file(args.pair_coeffs)

                equil_anneal_lmp.command("group solvate id <= {}".format(natoms_solvate_sys))
                equil_anneal_lmp.command("fix ic_prevention all momentum " +
                                         "{} linear 1 1 1 angular rescale".format(100))

                equil_anneal_lmp.command("dump dump_annealing all dcd {} {}".format(args.logsteps, equil_anneal_dcd))
                equil_anneal_lmp.command("dump_modify dump_annealing unwrap yes")
                equil_anneal_lmp.command("restart {} {} {}".format(args.logsteps*25,
                                                                   equil_anneal_rst,
                                                                   equil_anneal_rst))

                # thermo output
                equil_anneal_lmp.command("thermo_style custom " + " ".join(thermargs))
                equil_anneal_lmp.command("thermo_modify lost warn flush yes")
                equil_anneal_lmp.command("thermo {}".format(args.logsteps))

                # choose ensemble
                equil_anneal_lmp.command(("fix fix_nve_equil_anneal {} nve").format("all"))
                equil_anneal_lmp.command(("fix fix_temp_berendsen_equil_anneal {} temp/berendsen {} {} 0.1").format("all", start_anneal_temp, stop_anneal_temp))

                if args.lmps is not None:
                    equil_anneal_lmp.command(("fix fix_press_berendsen_equil_anneal {} press/berendsen iso 1.0 1.0 1.0").format("all"))

                # skip heating run if it is already finished but was restarted
                if rank == 0:
                    equil_anneal_llog = agul.LogUnification()
                    equil_anneal_llog.read_lmplog(equil_anneal_log)

                    # read last step or define it as None if there was None before
                    try:
                        last_step = equil_anneal_llog.data[-1]["Steps"][-1]
                    except IndexError:
                        last_step = None

                    del equil_anneal_llog

                last_step = comm.bcast(last_step, 0)

                if last_step < equil_steps:
                    equil_anneal_lmp.command("run {} upto".format(equil_steps))

                # short run at constant target temperature
                equil_anneal_lmp.command("unfix fix_temp_berendsen_equil_anneal")
                equil_anneal_lmp.command(("fix fix_temp_berendsen_equil_anneal {} temp/berendsen {} {} 0.1").format("all", stop_anneal_temp, stop_anneal_temp))

                if last_step < equil_steps+(equil_steps//5):
                    equil_anneal_lmp.command("run {}".format(equil_steps//5))

                del (equil_steps, last_step)
                equil_anneal_lmp.command("reset_timestep 0")

                #equil_anneal_lmp.command("unfix equilibration_ensemble_annealing")
                equil_anneal_lmp.command("unfix fix_nve_equil_anneal")
                equil_anneal_lmp.command("unfix fix_temp_berendsen_equil_anneal")

                if args.lmps is not None:
                    equil_anneal_lmp.command("unfix fix_press_berendsen_equil_anneal")

                equil_anneal_lmp.command("write_restart {}".format(equil_anneal_out))
                equil_anneal_lmp.command("clear")
                equil_anneal_lmp.close()

                if os.path.isfile(equil_anneal_rst) is True and rank == 0:
                    os.remove(equil_anneal_rst)

                #==============================#
                # equilibration: aggregate check
                #==============================#

                status_equil_anneal_agg = None

                if rank == 0:
                    solution_sys.import_dcd(equil_anneal_dcd)
                    solution_sys.read_frames(frame=-2, to_frame=-1)  # read last frame
                    solution_sys.close_dcd()

                    # refresh last coordinates and boxes
                    if args.lmps is not None:
                        solvate_sys.ts_coords.append(solution_sys.ts_coords[-1][:natoms_solvate_sys])
                        solvate_sys.ts_boxes.append(solution_sys.ts_boxes[-1])

                    status_equil_anneal_agg = check_aggregate(solvate_sys,
                                                              frame_id=-1,
                                                              atm_atm_dist=4,
                                                              debug=False)
                    #time.sleep(5)  # DEBUGGING

                status_equil_anneal_agg = comm.bcast(status_equil_anneal_agg, 0)

                if status_equil_anneal_agg is False:
                    # delete old folders
                    if rank == 0:
                        print("***Equil-Anneal Warning: Annealing failed! Starting all over.")
                        sl.rmtree(sysprep_dir)
                        sl.rmtree(quench_dir)
                        sl.rmtree(anneal_dir)

                    quench_success = False
                    anneal_success = False
                    time.sleep(5)  # DEBUGGING
                    continue

                del status_equil_anneal_agg

            #======================================#
            # annealing (solution) - productive run
            #======================================#

            anneal_steps = args.anneal_steps
            anneal_lmp = lammps()
            anneal_lmp.command("log {} append".format(anneal_dir+"anneal.lmplog"))

            if args.gpu:
                # neighbor list only on cpu or c_pe values are wrong
                anneal_lmp.command("package gpu 1 neigh no")
                anneal_lmp.command("suffix gpu")

            anneal_lmp.file(args.set)

            # no boundary conditions in vacuo
            if args.lmps is None:
                anneal_lmp.command("boundary f f f")

            # read restart from productive or from previous equilibration run
            if os.path.isfile(anneal_rst) is True:
                anneal_lmp.command("read_restart {}".format(anneal_rst))
            else:
                anneal_lmp.command("read_restart {}".format(equil_anneal_out))

            if args.pair_coeffs is not None:
                anneal_lmp.file(args.pair_coeffs)

            anneal_lmp.command("group solvate id <= {}".format(natoms_solvate_sys))
            anneal_lmp.command("fix ic_prevention all momentum " +
                               "{} linear 1 1 1 angular rescale".format(100))

            #anneal_lmp.command("dump dump_annealing all dcd {} {}".format(args.logsteps,
            #                                                              anneal_dcd))
            #anneal_lmp.command("dump_modify dump_annealing unwrap yes")
            anneal_lmp.command("restart {} {} {}".format(args.logsteps*25,
                                                         anneal_rst,
                                                         anneal_rst))

            # compute potential energy of the solvate (pair, bond, ...)
            # caveat: only possible with neigh no if calculating on the graphics
            #         card
            anneal_lmp.command("compute poteng_solvate all pe/atom")
            anneal_lmp.command("compute pe all reduce sum c_poteng_solvate")

            # extend thermo output
            anneal_lmp.command("thermo_style custom " + " ".join(thermargs) + " c_pe")
            anneal_lmp.command("thermo_modify lost warn flush yes")
            anneal_lmp.command("thermo {}".format(args.logsteps))

            if args.lmps is None or args.equil_anneal_ensemble == "nvt":
                anneal_lmp.command(("fix productive_ensemble_annealing {} "
                                    "nvt temp {} {} 0.1").format("all",
                                                                 stop_anneal_temp,
                                                                 stop_anneal_temp))
            else:
                anneal_lmp.command("group solvent id > {}".format(natoms_solvate_sys))
                anneal_lmp.command(("fix productive_ensemble_annealing {} "
                                    "npt temp {} {} 0.1 iso 1.0 1.0 1").format("all",
                                                                               stop_anneal_temp,
                                                                               stop_anneal_temp))

            #===============#
            # productive run
            #===============#
            total_anneal_runs = 5

            # get id of the last anneal run, all former dcd- and log-files
            if rank == 0:
                anneal_dcds = natsorted(["{}{}".format(anneal_dir, dcd_file) for dcd_file in
                                         os.listdir(anneal_dir) if
                                         re.match(r'\d+_anneal_\d+.dcd', dcd_file)])
                anneal_logs = natsorted(["{}{}".format(anneal_dir, log_file) for log_file in
                                         os.listdir(anneal_dir) if
                                         re.match(r'\d+_anneal_\d+.lmplog', log_file)])
                del (dcd_file, log_file)

                try:
                    # get last run number from file name
                    last_anneal_run = re.match(r'^\d+', os.path.basename(anneal_dcds[-1])).group(0)
                    # do not repeat the last run, continue with the next one
                    next_anneal_run = int(last_anneal_run) + 1
                    del last_anneal_run
                except IndexError:  # no file given (first start)
                    next_anneal_run = 0
            else:
                next_anneal_run = None

            #anneal_dcds = comm.bcast(anneal_dcds, 0)
            #anneal_logs = comm.bcast(anneal_logs, 0)
            next_anneal_run = comm.bcast(next_anneal_run, 0)

            if next_anneal_run == total_anneal_runs:
                next_anneal_run -= 1  # redo last aborted run

            for anneal_run in xrange(next_anneal_run, total_anneal_runs):
                # dcd stuff
                anneal_dcd = anneal_dir + "{}_anneal_{}".format(anneal_run, curcycle) + ".dcd"
                anneal_lmp.command("dump dump_annealing all dcd {} {}".format(args.logsteps, anneal_dcd))
                anneal_lmp.command("dump_modify dump_annealing unwrap yes")
                if anneal_dcd not in anneal_dcds:
                    anneal_dcds.append(anneal_dcd)

                # log stuff
                anneal_log = anneal_dir + "{}_anneal_{}".format(anneal_run, curcycle) + ".lmplog"
                anneal_lmp.command("log {}".format(anneal_log))
                if anneal_log not in anneal_logs:
                    anneal_logs.append(anneal_log)

                del (anneal_log)

                # perform current run
                anneal_lmp.command("run {}".format(anneal_steps))
                anneal_lmp.command("undump dump_annealing")
                anneal_lmp.command("write_restart {}".format(anneal_rst))

                # load last simulation frame
                status_anneal_agg = False
                normal_distributed_pe_quench = False
                pe = None

                if rank == 0:

                    #==========================================#
                    # productive: check aggregate of last frame
                    #==========================================#

                    solution_sys.import_dcd(anneal_dcd)
                    solution_sys.read_frames(frame=-2, to_frame=-1)  # read last frame
                    solution_sys.close_dcd()

                    # refresh last coordinates and boxes
                    if args.lmps is not None:
                        solvate_sys.ts_coords.append(solution_sys.ts_coords[-1][:natoms_solvate_sys])
                        solvate_sys.ts_boxes.append(solution_sys.ts_boxes[-1])

                    status_anneal_agg = check_aggregate(solvate_sys,
                                                        frame_id=-1,
                                                        atm_atm_dist=4,
                                                        debug=False)

                    #====================================================#
                    # productive: energy check of last 80 % of all frames
                    #====================================================#

                    energy_result = check_energy_convergence(anneal_logs,
                                                             keyword="c_pe",
                                                             perecentage=percentage_to_check,
                                                             debug=True)
                    skew               = energy_result[0]
                    p_skew             = energy_result[1]
                    p_normal           = energy_result[2]
                    #status_normal_dist = energy_result[3]
                    normal_distributed_pe_quench = energy_result[3]
                    min_pe             = energy_result[5]

                    # write results to the logfile
                    write_to_log("Skewness {} ".format(skew))
                    write_to_log("Pvalue (skewness) {} ".format(p_skew))
                    write_to_log("Pvalue (normality) {}\n\n".format(p_normal))
                    del (skew, p_skew, p_normal)

                del anneal_dcd
                status_anneal_agg = comm.bcast(status_anneal_agg, 0)
                normal_distributed_pe_quench = comm.bcast(normal_distributed_pe_quench, 0)

                #========================================#
                # check outcome of last annealing attempt
                #========================================#

                #DEBUGGING skip extended calculation
                anneal_success = True

                if status_anneal_agg is False:  # aggregate not ok
                    if rank == 0:
                        print("***Productive-Anneal Warning: Aggregate got dissolved! Annealing failed!")
                        time.sleep(5)
                    break
                elif normal_distributed_pe_quench is True:  # aggregate ok, energy ok
                    anneal_success = True
                    break
                else:  # aggregate ok, energy not yet
                    pass

            anneal_success = comm.bcast(anneal_success)
            anneal_lmp.command("clear")
            anneal_lmp.close()
            del (status_anneal_agg, anneal_steps)

            #=================================================================#
            # aftermath: extract frame with lowest pe of the solvate and check
            #            if its aggregate is intact
            #=================================================================#
            if anneal_success is True:
                #status_min_pe_anneal_agg = False

                if rank == 0:

                    # find dcd file with lowest pe/solvate
                    index_min_pe = None

                    # find minimum energy of all considered frames (i.e. all
                    # frames of an equilibrated system) and find the correspond-
                    # ing log file or rather its index which can be used for the
                    # dcd files to extract the coordinates of said frame
                    for index_logfile, anneal_log in enumerate(anneal_logs):
                        clog = agul.LogUnification()
                        clog.read_lmplog(anneal_log)

                        # get index of frame with lowest potential energy
                        if min_pe in clog.data[-1]["c_pe"]:
                            # last frame is current frame read
                            index_min_pe = clog.data[-1]["c_pe"].index(min_pe)
                            break

                    del (clog, min_pe)
                    solution_sys.ts_coords = []
                    solution_sys.import_dcd(anneal_dcds[index_logfile])
                    del index_logfile
                    solution_sys.read_frames(frame=index_min_pe,
                                             to_frame=index_min_pe + 1)
                    del index_min_pe
                    solution_sys.close_dcd()

                    #=============================#
                    # separate solute from solvent
                    #=============================#

                    if args.lmps is not None:
                        coords_solvent = solution_sys.ts_coords[-1][natoms_solvate_sys:]
                        solvent_sys.ts_coords.append(coords_solvent)
                        del coords_solvent
                        solvent_sys.write_xyz(solvent_anneal_out, -1)

                    # write anneal output
                    coords_solvate = solution_sys.ts_coords[-1][:natoms_solvate_sys]
                    solvate_sys.ts_coords.append(coords_solvate)

                    # final aggregate check (is everything o.k. with this frame?)
                    status_min_pe_agg = check_aggregate(solvate_sys, frame_id=-1, atm_atm_dist=4)

                    if status_min_pe_agg is True:
                        del coords_solvate
                        solvate_sys.write_xyz(solvate_anneal_out, -1)
                    else:
                        anneal_success = False

                status_min_pe_agg = comm.bcast(status_min_pe_agg, 0)
                anneal_success = comm.bcast(anneal_success, 0)
            else:
                # delete old folders
                if rank == 0:
                    print("***Anneal-Productive Info: Annealing finally failed! Starting all over.")
                    sl.rmtree(sysprep_dir)
                    sl.rmtree(quench_dir)
                    sl.rmtree(anneal_dir)

                quench_success = False
                continue

            #time.sleep(20)

        # tidy up
        if rank == 0:
            if args.lmps is not None:
                del solvent_sys
            del (solvate_sys, solution_sys)

    #==========================================================================#
    # 4. REQUENCHING
    #==========================================================================#
    if os.path.isfile(requench_dcd) is False:
        # solute
        solvate_sys = process_data(sysprep_out, quench_dcd)
        natoms_solvate_sys = len(solvate_sys.atoms)
        create_folder(requench_dir)

        #=========================================================#
        # read last frame and write temporary data file for lammps
        #=========================================================#
        if rank == 0:
            solvate_sys.read_xyz(solvate_anneal_out)
            solvate_sys.change_indices(incr=1, mode="increase")
            solvate_sys.write_lmpdat(tmp_solvate_anneal_out, -1,
                                     title="Input for requenching")
            solvate_sys.change_indices(incr=1, mode="decrease")

        # perform steepest descent minimization
        requench_lmp = lammps()
        requench_lmp.command("log {}".format(requench_log))
        requench_lmp.file(args.set)
        requench_lmp.command("read_data {}".format(tmp_solvate_anneal_out))

        if args.pair_coeffs is not None:
            requench_lmp.file(args.pair_coeffs)

        if rank == 0:
            os.remove(tmp_solvate_anneal_out)
            del tmp_solvate_anneal_out

        #requench_lmp.command("min_style sd")
        requench_lmp.command("min_style cg")
        #requench_lmp.command("fix minimize_requench all box/relax iso 1.0")
        requench_lmp.command("thermo_style custom " + " ".join(thermargs))
        requench_lmp.command("thermo_modify lost warn flush yes")
        requench_lmp.command("thermo {}".format(args.logsteps))
        requench_lmp.command("dump dump_requench all dcd {} {}".format(
            args.logsteps, requench_dcd))
        requench_lmp.command("dump_modify dump_requench unwrap yes")
        requench_lmp.command("neigh_modify every 1 delay 0 check yes")
        requench_steps = args.requench_steps
        requench_lmp.command("minimize 1e-6 1e-8 {} 100000".format(requench_steps))
        requench_lmp.command("clear")
        requench_lmp.close()

        # tidy up
        del requench_lmp
        if rank == 0:
            del solvate_sys

    if rank == 0:
        print("***Current cycle finished successfully!")

    time.sleep(2)
