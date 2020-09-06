import os
import numpy as np
import re
import glob
import shutil as sl
import argparse
from mpi4py import MPI
import ag_kwz as agk
import ag_lammps as aglmp
import ag_lammps_sim as aglmpsim
import time

# import ag_fileio
# import ag_lmplog as agl
# import ag_statistics as ags
import pdb

if __name__ == "__main__":
    # ==============================================================================#
    # Setup MPI
    # ==============================================================================#

    comm = MPI.COMM_WORLD
    size = comm.Get_size()  # number of processes in communicator
    rank = comm.Get_rank()  # process' id(s) within a communicator

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Place an aggregate into a solvent box and run a simulation",
    )

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
    solvent_pc_help = (
        "lammps'  script with lj, dreiding, etc. parameters for the solvent"
    )
    logsteps_help = (
        "log thermodynamic-, steps- and restart-files every" + "'logsteps' steps"
    )
    gpu_help = "Utilize lammps' GPU package (solvent relaxation not included)."
    solvent_gpu_help = "Use gpu for the solvent relaxation."
    cycles_help = "Number of aggregations."
    timeout_help = (
        "allowed duration of simulation, for resubmitting purposes;  should be < 24h"
    )
    pa_help = "The pattern in which looping order lmpa will be added, e.g. 0 1 2 3 3 1, repeats every 6 cycles"
    solution_paircoeffs_help = "Mixed pair types for the solvent and the solvate"
    dielectric_help = (
        "Reduce all coulombic interactions to 1/dielectric of their original strength."
    )

    # general
    parser.add_argument("-lmpm", default=None, metavar="*.lmpdat", help=lmpm_help)
    parser.add_argument(
        "-lmpa", default=None, nargs="*", metavar="*.lmpdat", help=lmpa_help
    )
    parser.add_argument("-lmps", default=None, metavar="*.lmpdat", help=lmps_help)
    parser.add_argument(
        "-lmps_dcd", default=None, metavar="*.lmpdat", help=lmps_dcd_help
    )
    # parser.add_argument("-solvate_resnames, metavar='cbz sac'")
    parser.add_argument("-set", metavar="*.lmpcfg", required=True, help=set_help)
    parser.add_argument(
        "-settings_solvent", metavar="*.lmpcfg", help=settings_solvent_help
    )
    parser.add_argument(
        "-pair_coeffs",
        default=None,
        metavar="*.lmpcfg",
        required=True,
        help=pair_coeffs_help,
    )
    parser.add_argument(
        "-solvent_paircoeffs", default=None, metavar="*.lmpcfg", help=solvent_pc_help
    )

    parser.add_argument(
        "-solution_paircoeffs",
        default=None,
        metavar="*.lmpcfg",
        help=solution_paircoeffs_help,
    )

    # parser.add_argument("-logsteps", type=int, default=1000, help=logsteps_help)
    parser.add_argument("-gpu", default=False, action="store_true", help=gpu_help)
    parser.add_argument(
        "-solvent_gpu", default=False, action="store_true", help=solvent_gpu_help
    )
    parser.add_argument("-cycles", type=int, default=5, help=cycles_help)
    parser.add_argument(
        "-timeout", metavar="00:01:00", default="00:00:05", help=timeout_help
    )
    parser.add_argument(
        "-pa", "-pattern_add", nargs="*", default=[0], type=int, help=pa_help
    )
    parser.add_argument("-dielectric", default=None, type=float, help=dielectric_help)

    args = parser.parse_args()

    if args.quench_np > size:
        raise Exception("More ranks selected than available")
