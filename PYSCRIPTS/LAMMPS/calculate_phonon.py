#!/usr/bin/env python

import os
from mpi4py import MPI
from lammps import lammps, PyLammps
from ag_lammps_sim import LmpSim

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()  # number of processes in communicator
RANK = COMM.Get_rank()  # process' id(s) within a communicator

def calculate_phonon(lmpcuts):
    lmp = lammps()
    # pylmp = PyLammps(ptr=lmp)
    # lmp.command("read_restart {}".format(lmpcuts.output_lmprst))
    lmp.file(lmpcuts.settings_file)
    lmp.command("box tilt large")
    lmp.command(f"read_data {lmpcuts.input_lmpdat}")
    lmpcuts.thermo(lmp)
    lmp.file(lmpcuts.pc_file)
    lmp.command("group phonon_IDs id 1:120")
    #phonon_cmd = "fix PHONON phonon_IDs phonon 20 1000 10000 GAMMA _ph"
    phonon_cmd = "fix PHONON all phonon 10 50000 1000000 GAMMA _ph nasr 50"
    lmp.command(phonon_cmd)
    lmpcuts.fix_hoover(lmp, "all", "nvt")
    lmpcuts.dump(lmp, unwrap=True)
    lmp.command("run {}".format(lmpcuts.runsteps))
    #lmp.command("dynamical_matrix all regular 0.000001")
    # use tools/phonon package for post-processing


if __name__ == "__main__":
    settings = "/hades/gadelmeier/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/md_settings/settings_dreiding_on.lmpcfg"
    paircoeffs = "/hades/gadelmeier/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/md_settings/CBZ_gaff-107_dreiding_on.lmpcfg"
    #q lmprst = "/hades/gadelmeier/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/3.unit_cell/1.cell_relax/CBZ_gaff-107/CBZIII/min_dreiding_on/CBZIII_gaff-107-min.lmprst"
    lmpdat = "/hades/gadelmeier/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/3.unit_cell/2.phonon/CBZIII/CBZIII_small_cell.lmpdat"
    #lmpdat = "/hades/gadelmeier/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/3.unit_cell/1.cell_relax/CBZ_gaff-107/CBZIII/CBZIII_gaff-107_supercell.lmpdat"

    lmpstuff = LmpSim(
        settings_file=settings,
        input_lmpdat=lmpdat,
        pc_file=paircoeffs,
        tstart=300,
        tstop=300,
        runsteps=1000000,
        logsteps=1000,
        output_dcd="out.dcd"
    )
    os.chdir(
        "/hades/gadelmeier/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/3.unit_cell/2.phonon/CBZIII/"
    )
    calculate_phonon(lmpstuff)
