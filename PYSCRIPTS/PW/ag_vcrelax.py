import copy
import os
import pathlib
import pdb
from typing import Iterable

import numpy as np

import ag_unify_md as agum


def find_pwout_w_lowest_energy(pathname, debug=False):
    """
    find_pwout_w_lowest_energy.
    Find pw output file which has the smallest energy of all files.

    Parameters
    ----------
    pathname : str
    Name of the path to search for pw output files.

    debug : bool
    Enable debugging mode.

    Returns
    -------
    file_w_lowest_energy : str
    Output file which has the lowest energy.
    """
    energy = 1e12

    for filename in pathlib.Path(".").glob("**/*.pwscf_out"):
        # get full path
        full_path = filename.resolve()
        pwscfout_file = full_path.as_posix()

        # read pwscf output
        pwsys = agum.Unification()
        pwsys.read_pwout(pwscfout_file)

        # get last energy from file (which should be the lowest)
        cenergy = pwsys.pw_other_info["ENERGIES"][-1]
        # print(pwscfout_file)
        # print(cenergy)

        # check if previous energy was higher
        if cenergy < energy:
            energy = cenergy
            file_w_lowest_energy = pwscfout_file

    return file_w_lowest_energy


def write_submit_script(
    pwin: str,
    pwout: str,
    slurm_commands: Iterable[str],
    pbs_commands: Iterable[str],
    scriptname: str = "default.bash",
) -> None:
    submit_script_part_1: str = """
    #  set variables according to cluster that is used
    if [ "$USER" == bccc34 ] || [ "$USER" == bccc013h ]; then
        module load intel64
        module load python/2.7-anaconda
        JOB_ID=$PBS_JOBID
        WORKING_DIR="$FASTTMP"
        SUBMIT_COMMAND="qsub"
        SUBMIT_DIR=${PBS_O_WORKDIR}
    elif [ "$USER" == gadelmeier ]; then
        JOB_ID=${SLURM_JOB_ID}
        WORKING_DIR="/scratch/"
        SUBMIT_DIR=$SLURM_SUBMIT_DIR
        SUBMIT_COMMAND="sbatch"
        export LD_LIBRARY_PATH="${HOME}/opt/compilers_and_libraries_2018/linux/mkl/lib/intel64/:${LD_LIBRARY_PATH}"
    fi
    """

    submit_script_part_2: str = """
    CUR_SCRIPT="$0"
    # directory where the pseudo potentials are listed
    export ESPRESSO_PSEUDO="$HOME/.qe_pps"
    """

    with open(scriptname, "w") as sn:
        # write the slurm commands first
        for key in slurm_commands:
            sn.write(f"#SBATCH {key}\n")

        for key in pbs_commands:
            sn.write(f"#PBS {key}\n")

        sn.write(submit_script_part_1)


write_submit_script("bla", "blabla")
