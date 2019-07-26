from __future__ import print_function, division
import os
import re
import glob
import shutil as sl
import argparse
from mpi4py import MPI
import ag_kwz as agk
import ag_lammps as aglmp
import ag_lammps_sim as aglmpsim
import time
#import ag_fileio
#import ag_lmplog as agl
#import ag_statistics as ags
import pdb

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
if __name__ == "__main__":
    percentage_to_check = 80
    #==============================================================================#
    # Argument Parsing
    #==============================================================================#
    parser = argparse.ArgumentParser(prog="kawska_zahn.py", formatter_class=argparse.RawTextHelpFormatter, description="Kawska-Zahn-Approach for accelerated crystallization simulations.")

    # arguments description
    lmpm_help = "Lammps' data-file of the main system."
    lmpa_help = """
    Lammps' data-file with atom-cube_sidetypes and coordinates for one single molecule to
    add to the current system. Atom types have to be defined already by the first
    data/restart file loaded!
    """
    lmps_help = "Create a solvent box for MD-Simulation. Data-file with one single solvent molecule to add."
    lmps_dcd_help = "DCD file of solvent box."
    set_help = "lammps' input-file/-script with basic simulation settings"
    settings_solvent_help = "lammps' input-file/-script with basic simulation settings for the solvent system"
    pair_coeffs_help = "lammps'  script with lj, dreiding, etc. parameters"
    solvent_pc_help = "lammps'  script with lj, dreiding, etc. parameters for the solvent"
    logsteps_help = "log thermodynamic-, steps- and restart-files every" + "'logsteps' steps"
    gpu_help = "utilize lammps' GPU package."
    cycles_help = "Number of aggregations."
    timeout_help = "allowed duration of simulation, for resubmitting purposes;  should be < 24h"
    pa_help = "The pattern in which looping order lmpa will be added, e.g. 0 1 2 3 3 1, repeats every 6 cycles"

    # general
    parser.add_argument("-lmpm", default=None, metavar="*.lmpdat", help=lmpm_help)
    parser.add_argument("-lmpa", default=None, nargs="*", metavar="*.lmpdat", help=lmpa_help)
    parser.add_argument("-lmps", default=None, metavar="*.lmpdat", help=lmps_help)
    parser.add_argument("-lmps_dcd", default=None, metavar="*.lmpdat", help=lmps_dcd_help)
    #parser.add_argument("-solvate_resnames, metavar='cbz sac'")
    parser.add_argument("-set", metavar="*.lmpcfg", required=True, help=set_help)
    parser.add_argument("-settings_solvent", metavar="*.lmpcfg", help=settings_solvent_help)
    parser.add_argument("-pair_coeffs", default=None, metavar="*.lmpcfg", required=True, help=pair_coeffs_help)
    parser.add_argument("-solvent_paircoeffs", default=None, metavar="*.lmpcfg", help=solvent_pc_help)
    parser.add_argument("-logsteps", type=int, default=1000, help=logsteps_help)
    parser.add_argument("-gpu", default=False, action="store_true", help=gpu_help)
    parser.add_argument("-cycles", type=int, default=5, help=cycles_help)
    parser.add_argument("-timeout", metavar="00:01:00", default="00:00:05", help=timeout_help)
    parser.add_argument("-pa", "-pattern_add", nargs="*", default=[0], type=int, help=pa_help)

    # quenching
    parser.add_argument("-quench_temp_start", type=int, default=5)
    parser.add_argument("-quench_temp_stop", type=int, default=5)
    parser.add_argument("-quench_steps", type=int, default=250000)
    parser.add_argument("-quench_logsteps", type=int, default=1000)

    # relax cut solvent
    parser.add_argument("-relax_cut_tstart", type=int, default=200)
    parser.add_argument("-relax_cut_tstop", type=int, default=250)
    parser.add_argument("-relax_cut_pstart", type=int, default=40)
    parser.add_argument("-relax_cut_pstop", type=int, default=10)
    parser.add_argument("-relax_cut_steps", type=int, default=50000)
    parser.add_argument("-relax_cut_logsteps", type=int, default=1000)

    # create voids in relaxed solvent
    parser.add_argument("-void_tstart", type=int, default=250)
    parser.add_argument("-void_tstop", type=int, default=300)
    parser.add_argument("-void_pstart", type=int, default=10)
    parser.add_argument("-void_pstop", type=int, default=1)
    parser.add_argument("-void_steps", type=int, default=2000)
    parser.add_argument("-void_logsteps", type=int, default=1000)

    # relax solvent in solution
    parser.add_argument("-relax_solv_tstart", type=int, default=250)
    #parser.add_argument("-relax_solv_style", type=str, default="berendsen")
    parser.add_argument("-relax_solv_tstop", type=int, default=300)
    parser.add_argument("-relax_solv_pstart", type=int, default=10)
    parser.add_argument("-relax_solv_pstop", type=int, default=1)
    parser.add_argument("-relax_solv_steps", type=int, default=2000)
    parser.add_argument("-relax_solv_logsteps", type=int, default=1000)

    # equilibrate solvate and solvent
    parser.add_argument("-heat_tstart", type=int, default=200)
    parser.add_argument("-heat_tstop", type=int, default=300)
    parser.add_argument("-heat_pstart", type=int, default=50)
    parser.add_argument("-heat_pstop", type=int, default=1)
    parser.add_argument("-heat_steps", type=int, default=50000)
    parser.add_argument("-heat_logsteps", type=int, default=1000)

    # annealing
    parser.add_argument("-anneal_tstart", type=int, default=300)
    parser.add_argument("-anneal_tstop", type=int, default=300)
    parser.add_argument("-anneal_pstart", type=int, default=1)
    parser.add_argument("-anneal_pstop", type=int, default=1)
    parser.add_argument("-anneal_steps", type=int, default=2000000)
    parser.add_argument("-anneal_steps_plus", type=int, default=500000)
    parser.add_argument("-anneal_logsteps", type=int, default=500)

    # requenching
    parser.add_argument("-requench_steps", type=int, default=150000)
    parser.add_argument("-requench_logsteps", type=int, default=500)

    args = parser.parse_args()

    #==============================================================================#
    # Remaining cycles and molecule to add pattern
    #==============================================================================#
    remaining_cycles, requench_out = agk.get_remaining_cycles(args.cycles)
    #resnames = set(args.solvate_resnames.split())

    #==============================================================================#
    # Kawska Zahn Approach
    #==============================================================================#
    PWD = os.getcwd()

    for curcycle in remaining_cycles:
        # get index of the pattern args.pa according to the current cycle
        pattern_idx = curcycle % len(args.pa)

        # get the index of args.lmpa according to the index the pattern gives
        # (args.pa is a list of indices (of args.lmpa) that are to be docked next)
        idx_lmpa = args.pa[pattern_idx]

        #==========================================================#
        # Define folders and files, retrieve stage of current cycle
        #==========================================================#
        if rank == 0:
            agk.write_to_log("Cycle: {:d}\n".format(curcycle))

        # declare folder names for each cycle
        sysprep_dir = PWD + "/sysprep_{}/".format(curcycle)
        quench_dir = PWD + "/quench_{}/".format(curcycle)
        anneal_dir = PWD + "/anneal_{}/".format(curcycle)
        requench_dir = PWD + "/requench_{}/".format(curcycle)

        # system preparation
        sysprep_out_lmpdat = sysprep_dir + "sysprep_out_{}.lmpdat".format(curcycle)
        lmpsettings_sysprep = aglmpsim.LmpSim(output_lmpdat=sysprep_out_lmpdat)

        # quench
        quench_out = quench_dir + "quench_out_{}.lmprst".format(curcycle)
        quench_rst = quench_dir + "quench_rst_{}.lmprst".format(curcycle)
        quench_dcd = quench_dir + "quench_{}.dcd".format(curcycle)
        quench_log = quench_dir + "quench_{}.lmplog".format(curcycle)
        lmpsettings_quench = aglmpsim.LmpSim(tstart=args.quench_temp_start, tstop=args.quench_temp_stop, logsteps=args.quench_logsteps, runsteps=args.quench_steps, pc_file=args.pair_coeffs, settings_file=args.set, input_lmpdat=sysprep_out_lmpdat, inter_lmprst=quench_rst, output_lmprst=quench_out, output_dcd=quench_dcd, output_lmplog=quench_log, gpu=args.gpu)

        # anneal -> solvent
        cut_solv_lmpdat = anneal_dir + "cut_solv_{}".format(curcycle) + "_out.lmpdat"
        cut_solv_rst = anneal_dir + "cut_solv_{}".format(curcycle) + "_tmp.rst"
        cut_solv_out = anneal_dir + "cut_solv_{}".format(curcycle) + "_out.lmprst"
        cut_solv_dcd = anneal_dir + "cut_solv_{}".format(curcycle) + ".dcd"
        cut_solv_log = anneal_dir + "cut_solv_{}".format(curcycle) + ".lmplog"
        lmpsettings_relax_cut = aglmpsim.LmpSim(tstart=args.relax_cut_tstart, tstop=args.relax_cut_tstop, pstart=args.relax_cut_pstart, pstop=args.relax_cut_pstop, logsteps=args.relax_cut_logsteps, runsteps=args.relax_cut_steps, pc_file=args.pair_coeffs, settings_file=args.settings_solvent, input_lmpdat=cut_solv_lmpdat, inter_lmprst=cut_solv_rst, output_lmprst=cut_solv_out, output_dcd=cut_solv_dcd, output_lmplog=cut_solv_log, gpu=args.gpu)

        void_solv_rst = anneal_dir + "void_solv_{}".format(curcycle) + "_tmp.rst"
        void_solv_out = anneal_dir + "void_solv_{}".format(curcycle) + "_out.lmprst"
        void_solv_dcd = anneal_dir + "void_solv_{}".format(curcycle) + ".dcd"
        void_solv_log = anneal_dir + "void_solv_{}".format(curcycle) + ".lmplog"
        lmpsettings_void = aglmpsim.LmpSim(tstart=args.void_tstart, tstop=args.void_tstop, pstart=args.void_pstart, pstop=args.void_pstop,logsteps=args.void_logsteps, runsteps=args.void_steps, pc_file=args.pair_coeffs,settings_file=args.settings_solvent, input_lmprst=lmpsettings_relax_cut.output_lmprst, inter_lmprst=void_solv_rst, output_lmprst=void_solv_out, output_dcd=void_solv_dcd,output_lmplog=void_solv_log, gpu=args.gpu)

        if args.lmps is not None:
            solution_lmpdat = anneal_dir + "solution_{}".format(curcycle) + "_out.lmpdat"
        else:
            solution_lmpdat = sysprep_out_lmpdat

        #relax_solv_in = anneal_dir + "relax_solv_{}".format(curcycle) + "_in.lmpdat"
        relax_solv_out = anneal_dir + "relax_solv_{}".format(curcycle) + "_out.lmprst"
        relax_solv_rst = anneal_dir + "relax_solv_{}".format(curcycle) + "_tmp.lmprst"
        relax_solv_dcd = anneal_dir + "relax_solv_{}".format(curcycle) + ".dcd"
        relax_solv_log = anneal_dir + "relax_solv_{}".format(curcycle) + ".lmplog"
        lmpsettings_relax_solv = aglmpsim.LmpSim(tstart=args.relax_solv_tstart, tstop=args.relax_solv_tstop, pstart=args.relax_solv_pstart, pstop=args.relax_solv_pstop,logsteps=args.relax_solv_logsteps, runsteps=args.relax_solv_steps, pc_file=args.pair_coeffs, settings_file=args.set,input_lmpdat=solution_lmpdat, inter_lmprst=relax_solv_rst,output_lmprst=relax_solv_out, output_dcd=relax_solv_dcd, output_lmplog=relax_solv_log,gpu=args.gpu)

        # anneal -> equilibration/heating
        heat_out = anneal_dir + "equil_anneal_{}".format(curcycle) + "_out.lmprst"
        heat_rst = anneal_dir + "equil_anneal_{}".format(curcycle) + "_tmp.lmprst"
        heat_dcd = anneal_dir + "equil_anneal_{}".format(curcycle) + ".dcd"
        heat_log = anneal_dir + "equil_anneal_{}".format(curcycle) + ".lmplog"
        lmpsettings_heat = aglmpsim.LmpSim(tstart=args.heat_tstart, tstop=args.heat_tstop, pstart=args.heat_pstart, pstop=args.heat_pstop, logsteps=args.heat_logsteps, runsteps=args.heat_steps,pc_file=args.pair_coeffs, settings_file=args.set, input_lmpdat=solution_lmpdat, input_lmprst=lmpsettings_relax_solv.output_lmprst, inter_lmprst=heat_rst, output_lmprst=heat_out, output_dcd=heat_dcd, output_lmplog=heat_log, gpu=args.gpu)

        if args.lmps is None:
            lmpsettings_heat.input_lmprst = lmpsettings_quench.output_lmprst
            #lmpsettings_heat.input_lmpdat = None
            #lmpsettings_heat.pstart = None
            #lmpsettings_heat.pstop = None

        # anneal -> productive
        anneal_out = anneal_dir + "anneal_{}".format(curcycle) + "_out.lmprst"
        anneal_rst = anneal_dir + "anneal_{}".format(curcycle) + "_tmp.lmprst"
        anneal_dcd = anneal_dir + "anneal_{}".format(curcycle) + ".dcd"
        anneal_log = anneal_dir + "anneal_{}".format(curcycle) + ".lmplog"
        lmpsettings_anneal = aglmpsim.LmpSim(tstart=args.anneal_tstart, tstop=args.anneal_tstop, pstart=args.anneal_pstart, pstop=args.anneal_pstop, logsteps=args.anneal_logsteps, runsteps=args.anneal_steps, pc_file=args.pair_coeffs, settings_file=args.set, input_lmpdat=solution_lmpdat, input_lmprst=lmpsettings_heat.output_lmprst, inter_lmprst=anneal_rst, output_lmprst=anneal_out, output_dcd=anneal_dcd, output_lmplog=anneal_log, gpu=args.gpu)

        # requench
        #tmp_solvate_anneal_out = requench_dir + "requench_{}".format(curcycle) + "_tmp_solvate_out.lmpdat"
        #requench_out = requench_dir + "requench_{}".format(curcycle) + "_out.lmpdat"
        requench_lmpdat = requench_dir + "requench_{}".format(curcycle) + ".lmpdat"
        requench_out = requench_dir + "requench_{}".format(curcycle) + "_out.lmprst"
        requench_rst = anneal_dir + "requench_{}".format(curcycle) + "_tmp.lmprst"
        requench_dcd = requench_dir + "requench_{}".format(curcycle) + ".dcd"
        requench_log = requench_dir + "requench_{}".format(curcycle) + ".lmplog"

        # not sure if "lmpsettings_requench = lmpsettings_anneal" is just a typo
        #lmpsettings_requench = lmpsettings_anneal = aglmpsim.LmpSim(logsteps=args.requench_logsteps, runsteps=args.requench_steps, pc_file=args.pair_coeffs, settings_file=args.set, input_lmpdat=solution_lmpdat, inter_lmprst=requench_rst, output_lmprst=requench_out, output_dcd=requench_dcd, output_lmplog=requench_log)
        lmpsettings_requench = aglmpsim.LmpSim(logsteps=args.requench_logsteps, runsteps=args.requench_steps, pc_file=args.pair_coeffs, settings_file=args.set, input_lmpdat=requench_lmpdat, inter_lmprst=requench_rst, output_lmprst=requench_out, output_dcd=requench_dcd, output_lmplog=requench_log)

        # important files from previous run
        pre_sysprep_out = "{0}/sysprep_{1}/sysprep_out_{1}.lmpdat".format(PWD, curcycle - 1)
        pre_solvent_anneal_out = "{0}/anneal_{1}/anneal_{0}_solvent_out.xyz".format(PWD, curcycle - 1)
        pre_requench_dcd = "{0}/requench_{1}/requench_{1}.dcd".format(PWD, curcycle - 1)

        quench_success = os.path.isfile(lmpsettings_quench.output_lmprst)
        anneal_success = os.path.isfile(lmpsettings_anneal.output_lmprst)

        #==========================================================================#
        # Aggregation
        #==========================================================================#

        # define main system
        if os.path.isfile(pre_requench_dcd) is True:
            main_prep_lmpdat = pre_sysprep_out
        else:
            main_prep_lmpdat = os.path.abspath(args.lmpm)

        while not os.path.isfile(lmpsettings_requench.output_lmprst):
            anneal_attempts = 0
            while not os.path.isfile(lmpsettings_anneal.output_lmprst):
                quench_attempts = 0
                while not os.path.isfile(lmpsettings_quench.output_lmprst):
                    sysprep_attempt = 0
                    while not os.path.isfile(lmpsettings_sysprep.output_lmpdat):
                        #==================================================================#
                        # 1. System Preparation
                        #==================================================================#
                        if rank == 0:
                            agk.create_folder(sysprep_dir)

                            # check if a previous run with a dcd file from requenching exists
                            if not os.path.isfile(pre_requench_dcd):
                                sysprep_success = agk.sysprep(lmpsettings_sysprep.output_lmpdat, main_prep_lmpdat, args.lmpa[idx_lmpa], dcd_add=None, frame_idx_main=-1)
                            else:
                                sysprep_success = agk.sysprep(lmpsettings_sysprep.output_lmpdat, main_prep_lmpdat, args.lmpa[idx_lmpa], dcd_main=pre_requench_dcd, dcd_add=None, frame_idx_main=-1)

                            if sysprep_success is False:
                                #agk.rename(sysprep_dir, sysprep_dir + "_failed_{}".format(sysprep_attempt))
                                agk.rename(sysprep_dir, sysprep_dir.rstrip("/") + "_failed")
                                sysprep_attempt += 1
                        else:
                            sysprep_success = None

                        sysprep_success = comm.bcast(sysprep_success, 0)

                        if sysprep_success is False and sysprep_attempt > 20:
                            exit(100)

                    #===================================================================
                    # 2. System Quenching
                    #===================================================================

                    if os.path.isfile(lmpsettings_quench.output_lmprst) is False:
                        if rank == 0:
                            agk.create_folder(quench_dir)

                        quench_success = agk.quench(lmpsettings_quench, main_prep_lmpdat)

                        if quench_success is False:
                            if rank == 0:
                                agk.rename(sysprep_dir, sysprep_dir.rstrip("/") + "_failed")
                                agk.rename(quench_dir, quench_dir.rstrip("/") + "_failed")

                            #del fail_appendix
                            quench_attempts += 1
                        else:
                            print("***Quenching-Info: Quenching done!")

                        # after 20 failed attempts, end run
                        if quench_attempts > 20 and quench_success is False:
                            exit(101)

                        del quench_attempts

                #======================================================================#
                # 3. ANNEALING
                #======================================================================#
                if os.path.isfile(lmpsettings_anneal.output_lmprst) is False:
                    if rank == 0:
                        agk.create_folder(anneal_dir)
                        solvate_sys = aglmp.read_lmpdat(lmpsettings_quench.input_lmpdat,
                                                        dcd=lmpsettings_quench.output_dcd)
                        solvate_sys_natoms = len(solvate_sys.atoms)
                        atm_idxs_solvate = range(solvate_sys_natoms)
                    else:
                        solvate_sys = None
                        atm_idxs_solvate = None

                    solvate_sys = comm.bcast(solvate_sys, 0)
                    atm_idxs_solvate = comm.bcast(atm_idxs_solvate, 0)

                    # check if solvent is needed
                    if args.lmps is not None:

                        # write data file for cut out solvent box
                        if rank == 0 and not os.path.isfile(cut_solv_lmpdat):
                            aglmp.cut_box(cut_solv_lmpdat, args.lmps, solvate_sys.ts_boxes[-1], args.lmps_dcd, frame_idx=-1)

                        # relax cut box
                        if not os.path.isdir(lmpsettings_relax_cut.output_lmprst):
                            lmpsettings_relax_cut.runsteps = 20000
                            # 20000 steps, 500 bar to 1 bar, 20 to 280 K, iso
                            agk.md_simulation(lmpsettings_relax_cut, group="all", style="berendsen", ensemble="nvt", keyword="iso")
                            lmpsettings_relax_cut.tstart = lmpsettings_relax_cut.tstop
                            agk.md_simulation(lmpsettings_relax_cut, group="all", style="berendsen", ensemble="npt", keyword="iso")

                        # create voids and write lammps data with solvate and solvent combined
                        if not os.path.isfile(lmpsettings_void.output_lmprst):
                            for _ in xrange(5):
                                no_clashes = agk.create_voids(lmpsettings_void, lmpsettings_quench.input_lmpdat, lmpsettings_quench.output_dcd)

                                if no_clashes is True:
                                    break

                            else:
                                print("Could not create large enough voids, void creation needs revision")
                                agk.rename(lmpsettings_void.output_lmprst, lmpsettings_void.inter_lmprst)
                                exit(102)

                        # combine solute and solvent
                        if rank == 0 and not os.path.isfile(solution_lmpdat):
                            aglmp.write_lmpdat(solution_lmpdat, sysprep_out_lmpdat, lmpdat_b=lmpsettings_relax_cut.input_lmpdat, dcd_a=lmpsettings_quench.output_dcd, dcd_b=lmpsettings_void.output_dcd, pair_coeffs=args.pair_coeffs)

                        # relax solvent molecules in solution, since solvent is always appended
                        # every atom id greater than the last one of the solvate has to be
                        # a solvent atom
                        agk.md_simulation(lmpsettings_relax_solv, group="group solvate id > {}".format(solvate_sys_natoms), style="berendsen", ensemble="nvt", keyword=None)

                    # ==========================================================
                    # heat the system to the  desired temperature

                    # carry out the complete run only if no 'output_lmprstt'
                    # file exists yet (intermediate restart files are used
                    # automatically by 'md_simulation' through LmpSim)
                    if not os.path.isfile(lmpsettings_heat.output_lmprst):
                        if args.lmps is not None:
                            agk.md_simulation(lmpsettings_heat, group="all", style="berendsen", ensemble="npt", keyword_min="iso", keyword="iso")
                        else:
                            # no solvent given -> keep volume constant
                            agk.md_simulation(lmpsettings_heat, group="all", style="berendsen", ensemble="nvt")

                    # check if aggregate is still ok
                    if rank == 0:
                        solution_sys = aglmp.read_lmpdat(lmpsettings_heat.input_lmpdat, lmpsettings_heat.output_dcd)
                        solution_sys_atoms_idxs = range(len(solution_sys.atoms))

                        # not sure if remains of previous code
                        #aggregate_ok = solution_sys.check_aggregate(solvate_sys, excluded_atm_idxs=solution_sys_atoms_idxs[solvate_sys_natoms:])

                        # check aggregation state of the last frame
                        aggregate_ok = solution_sys.check_aggregate(-1, excluded_atm_idxs=solution_sys_atoms_idxs[solvate_sys_natoms:])

                        # stop further calculations and start from the beginning
                        if not aggregate_ok:
                            if rank == 0:
                                agk.rename(sysprep_dir, sysprep_dir.rstrip("/") + "_failed")
                                agk.rename(quench_dir, quench_dir.rstrip("/") + "_failed")
                                agk.rename(anneal_dir, anneal_dir.rstrip("/") + "_failed")

                    else:
                        aggregate_ok = False

                    aggregate_ok = comm.bcast(aggregate_ok, 0)

                    # reset all failed attempts and restart
                    if not aggregate_ok:
                        anneal_attempts = 0
                        quench_attempts = 0
                        sysprep_attempt = 0
                        continue

                    # productive run
                    if args.lmps is not None:
                        anneal_success, anneal_dcds, log_files = agk.anneal_productive(lmpsettings_anneal, atm_idxs_solvate, percentage_to_check, "npt", keyword="iso")
                    else:
                        anneal_success, anneal_dcds, log_files = agk.anneal_productive(lmpsettings_anneal, atm_idxs_solvate, percentage_to_check, "nvt", output=anneal_dir + "anneal_{}".format(curcycle))

                    # start all over if annealing failed
                    if not anneal_success:
                        anneal_attempts = 0
                        quench_attempts = 0
                        sysprep_attempt = 0
                        agk.rename(sysprep_dir, sysprep_dir.rstrip("/") + "_failed")
                        agk.rename(quench_dir, quench_dir.rstrip("/") + "_failed")
                        agk.rename(anneal_dir, anneal_dir.rstrip("/") + "_failed")
                        #pdb.set_trace()
                        continue

            #======================================================================#
            # 4. REQUENCHING
            #======================================================================#
            if rank == 0:
                agk.create_folder(requench_dir)

                # write input lmpdat with the best frame from the annealing
                if os.path.isfile(lmpsettings_requench.input_lmpdat) is False:
                    # gather all lmplog and dcd files
                    anneal_lmplog_files = glob.glob(r"{}/*[0-9]_anneal_[0-9]*.lmplog".format(anneal_dir))
                    anneal_dcds = glob.glob(r"{}/*[0-9]_anneal_[0-9]*.dcd".format(anneal_dir))
                    # find frame which scores best (dcd and idx) and its value
                    best_dcd, best_idx, best_val = agk.find_best_frame(anneal_lmplog_files, anneal_dcds, thermo="c_pe_solvate_complete", percentage_to_check=percentage_to_check)
                    # write the data file for requenching
                    agk.write_requench_data(lmpsettings_sysprep.output_lmpdat, best_dcd, best_idx, output_lmpdat_a=lmpsettings_requench.input_lmpdat)

            requench_success = agk.requench(lmpsettings_requench)

            if requench_success is False:
                anneal_attempts = 0
                quench_attempts = 0
                sysprep_attempt = 0
                agk.rename(sysprep_dir, sysprep_dir.rstrip("/") + "_failed")
                agk.rename(quench_dir, quench_dir.rstrip("/") + "_failed")
                agk.rename(anneal_dir, anneal_dir.rstrip("/") + "_failed")
                agk.rename(requench_dir, requench_dir.rstrip("/") + "_failed")
                continue

            print("***Requenching-Info: Requenching done!")
