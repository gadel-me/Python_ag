from __future__ import print_function, division
import os
import re
import shutil as sl
import argparse
from mpi4py import MPI
import ag_kwz as agk
import ag_statistics as ags
import vmd

"""
Kawska-Zahn approach to aggregate crystals.

This script is doing a Kawska-Zahn approach to crystallize given molecules using
lammps as the driver for molecular dynamics simulations. Equilibration is
checked and the simulation time is elongated if the system has not equilibrated
yet.

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

#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator

#==============================================================================#
# Helper functions and variables
#==============================================================================#
PWD = os.getcwd()


def get_finished_cycles():
    """
    Gather all finished cycles.

    Returns
    -------
    finished_cycles : list of ints
        all ids of the cycles that were finished
    """
    finished_cycles = []
    folders_pwd = ["{}/{}".format(PWD, i) for i in os.listdir(PWD) if os.path.isdir(i)]

    # get last cycle from directory
    for folder in folders_pwd:
        cycle = re.match(r'.*?([0-9]+)$', folder).group(1)
        cycle = int(cycle)

        # avoid duplicates
        if cycle not in finished_cycles:
            finished_cycles.append(cycle)

    return finished_cycles


def get_next_cycle(finished_cycles):
    """
    requench_out does not exist -> current cycle has not finished yet
    """
    if len(finished_cycles) >= 1:
        current_cycle = max(finished_cycles)

        if os.path.isfile(requench_out) is True:
            next_cycle = current_cycle + 1
        else:
            next_cycle = current_cycle

        #del requench_out
    else:
        next_cycle = 0

    return (current_cycle, next_cycle)


def get_remaining_cycles():
    """
    Scan current folder names for numbers to get the current cycle.

    Returns
    -------
    remaining_cycles : int
        remaining cycles

    """
    finished_cycles = get_finished_cycles()
    current_cycle, next_cycle = get_next_cycle(finished_cycles)
    requench_out = PWD + "/requench_{}/".format(current_cycle) + "requench_{}.dcd".format(current_cycle)
    total_cycles = range(args.cycles)
    idx_next_cycle = total_cycles.index(next_cycle)
    del total_cycles

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
    return (remaining_cycles, requench_out)


def create_folder(folder):
    """
    Create folder or skip creation if it already exists
    """
    try:
        os.mkdir(folder, 0755)
    except OSError:
        print("***Info: Folder {} already exists!".format(folder))


def write_to_log(string, filename="kwz_log"):
    """
    Write string 's' to log file named 'kwz.log'.
    """
    with open("kwz_log", "a") as kwz_log:
        kwz_log.write(string)

if __name__ == "__main__":
    #==============================================================================#
    # Argument Parsing
    #==============================================================================#
    parser = argparse.ArgumentParser(prog="kawska_zahn.py", formatter_class=argparse.RawTextHelpFormatter, description="Kawska-Zahn-Approach for accelerated crystallization simulations.")

    # general
    parser.add_argument("-lmpm", default=None, metavar="*.lmpdat", help="Lammps' data-file of the main system.")
    parser.add_argument("-lmpa", default=None, nargs="*", metavar="*.lmpdat", help="Lammps' data-file with atom-cube_sidetypes and coordinates for one single molecule to add to the current system. Atom types have to be defined already by the first data/restart file loaded!")
    parser.add_argument("-lmps", default=None, metavar="*.lmpdat", help="Create a solvent box for MD-Simulation. Data-file with one single solvent molecule to add.")
    parser.add_argument("-lmps_dcd", default=None, metavar="*.lmpdat", help="DCD file of solvent box.")
    parser.add_argument("-set", metavar="*.lmpcfg", required=True, help="lammps' input-file/-script with basic simulation " + "settings")
    parser.add_argument("-settings_solvent", metavar="*.lmpcfg", help="lammps' input-file/-script with basic simulation " +      "settings for the solvent system")
    parser.add_argument("-pair_coeffs", default=None, metavar="*.lmpcfg", required=True, help="lammps'  script with lj, dreiding, etc. parameters")
    parser.add_argument("-solvent_paircoeffs", default=None, metavar="*.lmpcfg", help="lammps'  script with lj, dreiding, etc. parameters " + "for the solvent")
    parser.add_argument("-logsteps", type=int, default=1000, help="log thermodynamic-, steps- and restart-files every" + "'logsteps' steps")
    parser.add_argument("-gpu", default=False, action="store_true", help="utilize lammps' GPU package.",)
    parser.add_argument("-cycles", type=int, default=5, help="Number of aggregations.")
    parser.add_argument("-timeout", metavar="00:01:00", default="00:00:05", help="allowed duration of simulation, for resubmitting purposes;  should be < 24h")
    parser.add_argument("-pa", "-pattern_add", nargs="*", default=0, type=int, help="The pattern in which looping order lmpa will be added, e.g. 0 1 2 3 3 1, repeats every 6 cycles")

    # quenching
    parser.add_argument("-quench_temp_start", type=int, default=5)
    parser.add_argument("-quench_temp_stop", type=int, default=5)
    parser.add_argument("-quench_steps", type=int, default=250000)
    parser.add_argument("-quench_logsteps", type=int, default=1000)

    # cut and relax solvent
    parser.add_argument("-cut_tstart", type=int, default=200)
    parser.add_argument("-cut_tstop", type=int, default=250)
    parser.add_argument("-cut_pstart", type=int, default=40)
    parser.add_argument("-cut_pstop", type=int, default=10)
    parser.add_argument("-cut_steps", type=int, default=50000)

    # create voids in relaxed solvent
    parser.add_argument("-void_tstart", type=int, default=250)
    parser.add_argument("-void_tstop", type=int, default=300)
    parser.add_argument("-void_pstart", type=int, default=10)
    parser.add_argument("-void_pstop", type=int, default=1)
    parser.add_argument("-void_steps", type=int, default=2000)
    parser.add_argument("-void_logsteps", type=int, default=1000)

    # equilibrate solvate and solvent
    parser.add_argument("-start_relax_solvent_temp", type=int, default=10)
    parser.add_argument("-stop_relax_solvent_temp", type=int, default=10)
    parser.add_argument("-relax_solvent_steps", type=int, default=50000)
    parser.add_argument("-start_equil_anneal_temp", type=int, default=10)
    parser.add_argument("-stop_equil_anneal_temp", type=int, default=300, help="Temperature at which system shall be annealed")
    parser.add_argument("-equil_anneal_steps", type=int, default=500000)
    parser.add_argument("-equil_anneal_ensemble", type=str, default="npt", help="nvt or npt")
    parser.add_argument("-anneal_steps", type=int, default=2000000)
    parser.add_argument("-additional_anneal_steps", type=int, default=100000)
    parser.add_argument("-anneal_logsteps", type=int, default=500)
    parser.add_argument("-requench_steps", type=int, default=150000)

    args = parser.parse_args()

    #==============================================================================#
    # Remaining cycles and molecule to add pattern
    #==============================================================================#
    remaining_cycles, requench_out = get_remaining_cycles()


    #==============================================================================#
    # Kawska Zahn Approach
    #==============================================================================#

    for curcycle, idx_lmpa in remaining_cycles:

        #==========================================================#
        # Define folders and files, retrieve stage of current cycle
        #==========================================================#
        if rank == 0:
            write_to_log("Cycle: {:d}\n".format(curcycle))

        # declare folder names for each cycle
        sysprep_dir = PWD + "/sysprep_{}/".format(curcycle)
        quench_dir = PWD + "/quench_{}/".format(curcycle)
        anneal_dir = PWD + "/anneal_{}/".format(curcycle)
        requench_dir = PWD + "/requench_{}/".format(curcycle)

        # system preparation
        sysprep_out = sysprep_dir + "sysprep_out_{}.lmpdat".format(curcycle)

        # quench
        quench_out = quench_dir + "quench_out_{}.lmprst".format(curcycle)
        quench_rst = quench_dir + "quench_rst_{}.lmprst".format(curcycle)
        quench_dcd = quench_dir + "quench_{}.dcd".format(curcycle)
        quench_log = quench_dir + "quench_{}.lmplog".format(curcycle)
        lmpcuts_quench = agk.LmpShortcuts(tstart=args.quench_temp_start, tstop=args.quench_temp_stop, logsteps=args.quench_logsteps, runsteps=args.quench_steps, pc_file=args.pair_coeffs, settings_file=args.set, input_lmpdat=sysprep_out, inter_lmprst=quench_rst, output_lmprst=quench_out, output_dcd=quench_dcd, output_lmplog=quench_log, gpu=args.gpu)

        # anneal -> solvent
        cut_solv = anneal_dir + "cut_solv_{}".format(curcycle) + "_out.lmpdat"
        cut_solv_rst = anneal_dir + "cut_solv_{}".format(curcycle) + "_tmp.rst"
        cut_solv_out = anneal_dir + "cut_solv_{}".format(curcycle) + "_out.lmprst"
        cut_solv_dcd = anneal_dir + "cut_solv_{}".format(curcycle) + ".dcd"
        cut_solv_log = anneal_dir + "cut_solv_{}".format(curcycle) + ".lmplog"
        lmpcuts_cut_solv = agk.LmpShortcuts(tstart=args.cut_tstart, tstop=args.cut_tstop, pstart=args.cut_pstart, pstop=args.cut_pstop, logsteps=args.void_logsteps, runsteps=args.void_steps, pc_file=args.pair_coeffs, settings_file=args.settings_solvent, input_lmpdat=cut_solv, inter_lmprst=cut_solv_rst, output_lmprst=cut_solv_out, output_dcd=cut_solv_dcd, output_lmplog=cut_solv_log, gpu=args.gpu)

        void_solvate_solvent = anneal_dir + "void_solv_{}".format(curcycle) + "_out.lmpdat"
        void_solv_rst = anneal_dir + "void_solv_{}".format(curcycle) + "_tmp.rst"
        void_solv_out = anneal_dir + "void_solv_{}".format(curcycle) + "_out.lmprst"
        void_solv_dcd = anneal_dir + "void_solv_{}".format(curcycle) + ".dcd"
        void_solv_log = anneal_dir + "void_solv_{}".format(curcycle) + ".lmplog"
        lmpcuts_void = agk.LmpShortcuts(tstart=args.void_tstart, tstop=args.void_tstop, pstart=args.void_pstart, pstop=args.void_pstop, logsteps=args.void_logsteps, runsteps=args.void_steps, pc_file=args.pair_coeffs, settings_file=args.settings_solvent, input_lmpdat=cut_solv_out, inter_lmprst=void_solv_rst, output_lmpdat=void_solvate_solvent, output_lmprst=void_solv_out, output_dcd=void_solv_dcd, output_lmplog=void_solv_log, gpu=args.gpu)

        relax_solv_in = anneal_dir + "relax_solv_{}".format(curcycle) + "_in.lmpdat"
        relax_solv_out = anneal_dir + "relax_solv_{}".format(curcycle) + "_out.lmprst"
        relax_solv_rst = anneal_dir + "relax_solv_{}".format(curcycle) + "_tmp.lmprst"
        relax_solv_dcd = anneal_dir + "relax_solv_{}".format(curcycle) + ".dcd"
        relax_solv_log = anneal_dir + "relax_solv_{}".format(curcycle) + ".lmplog"
        lmpcuts_relax = agk.LmpShortcuts(tstart=args.void_tstart, tstop=args.void_tstop, logsteps=args.void_logsteps, runsteps=args.void_steps, pc_file=args.pair_coeffs, settings_file=args.settings_solvent, input_lmpdat=cut_solv_out, inter_lmprst=void_solv_rst, output_lmpdat=void_solvate_solvent, output_lmprst=void_solv_out, output_dcd=void_solv_dcd, output_lmplog=void_solv_log, gpu=args.gpu)

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
        #requench_out          = requench_dir + "requench_{}".format(curcycle) + "_out.lmpdat"
        requench_dcd           = requench_dir + "requench_{}".format(curcycle) + ".dcd"
        requench_log           = requench_dir + "requench_{}".format(curcycle) + ".lmplog"

        # important files from previous run
        pre_sysprep_out = "{0}/sysprep_{1}/sysprep_out_{1}.lmpdat".format(PWD, curcycle - 1)
        pre_solvent_anneal_out = "{0}/anneal_{1}/anneal_{0}_solvent_out.xyz".format(PWD, curcycle - 1)
        pre_requench_dcd = "{0}/requench_{1}/requench_{1}.dcd".format(PWD, curcycle - 1)

        quench_success = os.path.isfile(quench_out)
        anneal_success = os.path.isfile(solvate_anneal_out)

        #==========================================================================#
        # Aggregation
        #==========================================================================#

        # define main system
        if os.path.isfile(pre_requench_dcd) is True:
            main_prep_lmpdat = pre_sysprep_out
        else:
            main_prep_lmpdat = os.path.abspath(args.lmpm)

        while anneal_success is False:
            anneal_attempts = 0
            while os.path.isfile(quench_out) is False:
                quench_attempts = 0
                while os.path.isfile(sysprep_out) is False:
                    sysprep_attempts = 0
                    #==================================================================#
                    # 1. System Preparation
                    #==================================================================#
                    if rank == 0:
                        create_folder(sysprep_dir)
                        sysprep_success = agk.sysprep(sysprep_out, main_prep_lmpdat, args.lmpa[idx_lmpa], dcd_add=requench_dcd, frame_idx=-1)

                        if sysprep_success is False:
                            sl.move(sysprep_dir, sysprep_dir + "failed_{}".format(sysprep_attempts))
                            sysprep_attempts += 1

                    if sysprep_success is False and sysprep_attempts > 20:
                        exit(100)

                #===================================================================
                # 2. System Quenching
                #===================================================================

                if os.path.isfile(quench_out) is False:
                    if rank == 0:
                        create_folder(quench_dir)
                    quench_success = agk.quench(lmpcuts_quench, main_prep_lmpdat)

                    if quench_success is False:
                        fail_appendix = "failed_{}".format(quench_attempts)
                        sl.move(sysprep_dir, sysprep_dir.rstrip("/") + fail_appendix)
                        sl.move(quench_dir, quench_dir.rstrip("/") + fail_appendix)
                        del fail_appendix
                        quench_attempts += 1
                    else:
                        print("***Quenching-Info: Quenching done!")
                        del quench_attempts

                    if quench_attempts > 20 and quench_success is False:
                        exit(101)

            #======================================================================#
            # 3. ANNEALING
            #======================================================================#
            if os.path.isfile(solvate_anneal_out) is False:
                if rank == 0:
                    create_folder(anneal_dir)

                solvate_sys = agk.read_system(sysprep_out, dcd=lmpcuts_quench.output_dcd)

                # check if solvent is needed
                if args.lmps is not None:

                    # write data file for cut out solvent box
                    if not os.path.isfile(cut_solv):
                        if rank == 0:
                            agk.cut_box(cut_solv_out, args.lmps, solvate_sys.ts_boxes[-1], args.lmps_dcd, frame_idx=-1)

                    if not os.path.isdir(cut_solv_out):
                        agk.relax_box(lmpcuts_cut_solv)

                    # create voids and write lammps data with solvate and solvent combined
                    if not os.path.isfile(void_solvate_solvent):
                        agk.create_voids(lmpcuts_void, sysprep_out, lmpcuts_quench.output_dcd)

                    if not os.path.isfile(relax_solv_out):
                        agk.relax_box(lmpcuts_relax)
