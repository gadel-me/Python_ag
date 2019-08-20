#!/usr/bin/env python
from __future__ import print_function, division
import pdb
import os
import argparse
import numpy as np
from mpi4py import MPI
from lammps import lammps, PyLammps
#import ag_lammps as aglmp
import ag_lammps_sim as aglmpsim

"""
This script intends to make a run for the so-called 'sandwich' system, i.e. form II
of carbamazepine with solvent inside its cavities and solvent above and below
as well. This calculation will test if water and thf stay inside the 'nano tubes'
or if they will leave them over time.
"""


#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator


def relax_group(lmpcuts, ensemble, group="all", keyword=None, create_velocity=True):
    """
    Relax the group of atoms of the given system.

    Parameters
    ----------
    lmpcuts : ag_lammps_sim.LmpSim
        basic settings
    ensemble : str ("npt" or "nvt")
        the type of simulation
    keyword : str or None
        one of lammps' keywords for npt runs

    """
    lmp = lammps()
    pylmp = PyLammps(ptr=lmp)
    lmp.command("log {} append".format(lmpcuts.output_lmplog))
    #lmpcuts.use_gpu(lmp, neigh=False)

    lmp.file(lmpcuts.settings_file)
    lmp.command("box tilt large")  # ignore too tilted boxes

    lmpcuts.load_system(lmp)
    lmpcuts.dump(lmp)
    lmpcuts.thermo(lmp)

    if lmpcuts.pc_file is not None:
        lmp.file(lmpcuts.pc_file)

    # define the atoms that may move during the simulation
    if group == "all":
        lmpcuts.fix_hoover(lmp, group, ensemble, keyword)
    else:
        lmp.command("group group_1 {}".format(group))
        lmpcuts.fix_hoover(lmp, "group_1", ensemble, keyword)

    lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")

    # pre-optimization
    if lmpcuts.input_lmprst is None:
        lmp.command("min_style cg")
        lmp.command("min_modify dmax 0.5")
        lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

    if create_velocity is True:
        lmp.command("velocity all create {} {} mom yes rot yes dist gaussian".format(lmpcuts.tstart, np.random.randint(29847587)))

    lmp.command("run {}".format(lmpcuts.runsteps))

    # tidy up before closing
    lmpcuts.unfix_undump(pylmp, lmp)
    #lmp.command("reset_timestep 0")
    lmp.command("write_restart {}".format(lmpcuts.output_lmprst))
    lmp.command("clear")

    # close lammps
    lmp.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-lmpdat", default=None)
    #parser.add_argument("-lmpdat_CBZII_only", default=None, help="Data file with solely CBZ")
    parser.add_argument("-relax_thread", action="store_true", default=False)
    parser.add_argument("-cv", action="store_true", default=True, help="Create an initial velocity (i.e. neglect existing one)")

    # relax solvent
    parser.add_argument("-solvent_temp_start", type=int, default=295)
    parser.add_argument("-solvent_temp_stop", type=int, default=295)
    parser.add_argument("-solvent_steps", type=int, default=100000)
    parser.add_argument("-solvent_logsteps", type=int, default=1000)
    parser.add_argument("-solvent_pstart", type=int, default=50)
    parser.add_argument("-solvent_pstop", type=int, default=1)
    parser.add_argument("-solvent_group", default="all")

    # relax thread
    parser.add_argument("-thread_temp_start", type=int, default=295)
    parser.add_argument("-thread_temp_stop", type=int, default=295)
    parser.add_argument("-thread_steps", type=int, default=500000)
    parser.add_argument("-thread_logsteps", type=int, default=1000)
    parser.add_argument("-thread_group", default="all")

    # relax all
    parser.add_argument("-all_temp_start", type=int, default=295)
    parser.add_argument("-all_temp_stop", type=int, default=295)
    parser.add_argument("-all_steps", type=int, default=1000000)
    parser.add_argument("-all_logsteps", type=int, default=1000)
    parser.add_argument("-all_pstart", type=int, default=1)
    parser.add_argument("-all_pstop", type=int, default=1)
    parser.add_argument("-all_group", default="all")

    # general settings
    parser.add_argument("-set", metavar="*.lmpcfg", required=True)
    parser.add_argument("-pair_coeffs", default=None, metavar="*.lmpcfg")
    parser.add_argument("-timeout", metavar="00:01:00", default="00:00:05")

    args = parser.parse_args()

    # relax solvent
    lmpsetting_relax_solvent = aglmpsim.LmpSim(
        tstart=args.solvent_temp_start,
        tstop=args.solvent_temp_stop,
        pstart=args.solvent_pstart,
        pstop=args.solvent_pstop,
        logsteps=args.solvent_logsteps,
        runsteps=args.solvent_steps,
        pc_file=args.pair_coeffs,
        settings_file=args.set,
        input_lmpdat=args.lmpdat,
        inter_lmprst="CBZII_relax_solvent_inter.lmprst",
        output_lmprst="CBZII_relax_solvent_out.lmprst",
        output_dcd="CBZII_relax_solvent.dcd",
        output_lmplog="CBZII_relax_solvent.lmplog")

    # relax thread
    lmpsetting_relax_thread = aglmpsim.LmpSim(
        tstart=args.thread_temp_start,
        tstop=args.thread_temp_stop,
        logsteps=args.thread_logsteps,
        runsteps=args.thread_steps,
        pc_file=args.pair_coeffs,
        settings_file=args.set,
        input_lmprst=lmpsetting_relax_solvent.output_lmprst,
        inter_lmprst="CBZII_relax_thread_inter.lmprst",
        output_lmprst="CBZII_relax_thread_out.lmprst",
        output_dcd="CBZII_relax_thread.dcd",
        output_lmplog="CBZII_relax_thread.lmplog")

    # relax all
    lmpsetting_relax_all = aglmpsim.LmpSim(
        tstart=args.all_temp_start,
        tstop=args.all_temp_stop,
        pstart=args.all_pstart,
        pstop=args.all_pstop,
        logsteps=args.all_logsteps,
        runsteps=args.all_steps,
        pc_file=args.pair_coeffs,
        settings_file=args.set,
        input_lmprst=lmpsetting_relax_thread.output_lmprst,
        inter_lmprst="CBZII_relax_all_inter.lmprst",
        output_lmprst="CBZII_relax_all_out.lmprst",
        output_dcd="CBZII_relax_all.dcd",
        output_lmplog="CBZII_relax_all.lmplog")

    # relax the solvent first
    if not os.path.isfile(lmpsetting_relax_solvent.output_lmprst):
        relax_group(lmpsetting_relax_solvent, "nvt", group=args.solvent_group, create_velocity=args.cv)

    # relax the thread in an nvt environment second
    if args.relax_thread is True and not os.path.isfile(lmpsetting_relax_all.output_lmprst):
        relax_group(lmpsetting_relax_thread, "nvt", group=args.thread_group, create_velocity=args.cv)

    # relax the whole system
    if not os.path.isfile(lmpsetting_relax_all.output_lmprst):

        if args.relax_thread is False:
            lmpsetting_relax_all.input_lmprst = lmpsetting_relax_solvent.output_lmprst
        #else:
        #    lmpsetting_relax_all.input_lmprst = lmpsetting_relax_thread.output_lmprst

        relax_group(lmpsetting_relax_all, "npt", group=args.all_group, keyword="tri", create_velocity=args.cv)
