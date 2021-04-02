from simulate_sandwich import relax_group

import argparse
import numpy as np
from mpi4py import MPI
from lammps import lammps, PyLammps

# import ag_lammps as aglmp
import ag_lammps_sim as aglmpsim

"""
This script intends to make a run for the so-called 'sandwich' system, i.e. form II
of carbamazepine with solvent inside its cavities and solvent above and below
as well. This calculation will test if water and thf stay inside the 'nano tubes'
or if they will leave them over time.
"""


# ==============================================================================#
# Setup MPI
# ==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator

lmpsettings = aglmpsim.LmpSim(
    tstart=300,
    tstop=300,
    pstart=1.0,
    pstop=1.0,
    logsteps=1000,
    runsteps=1000000,
    pc_file="/hades/gadelmeier/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/md_settings/CBZ_gaff-107_dreiding_on.lmpcfg",
    settings_file="/hades/gadelmeier/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/md_settings/standard_settings_dreiding_on.lmpcfg",
    input_lmpdat="/hades/gadelmeier/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/3.unit_cell/1.cell_relax/CBZ_gaff-107/CBZII/npt_2000K-2000K_dreiding_on-1/CBZII_gaff-107-npt_2000K-2000K-1.lmpdat",
    inter_lmprst="CBZII_107_dreiding_on_npt_300K-300K_dreiding_on_berendsen-melt.lmprst",
    output_lmprst="CBZII_107_dreiding_on_npt_300K-300K_dreiding_on_berendsen-melt.lmprst",
    output_dcd="CBZII_107_dreiding_on_npt_300K-300K_dreiding_on_berendsen-melt.dcd",
    output_lmplog="CBZII_107_dreiding_on_npt_300K-300K_dreiding_on_berendsen-melt.lmplog",
)

lmp = lammps()
pylmp = PyLammps(ptr=lmp)
lmp.command("log {} append".format(lmpsettings.output_lmplog))
# lmpsettings.use_gpu(lmp, neigh=False)

lmp.file(lmpsettings.settings_file)
lmp.command("box tilt large")  # ignore too tilted boxes

lmpsettings.load_system(lmp)
lmpsettings.dump(lmp)
lmpsettings.thermo(lmp)

if lmpsettings.pc_file is not None:
    lmp.file(lmpsettings.pc_file)

lmpsettings.fix_berendsen(lmp, "all", "npt", "aniso")
lmp.command("fix ic_prevention all momentum 100 linear 1 1 1 angular rescale")

# pre-optimization
if lmpsettings.input_lmprst is None:
    lmp.command("min_style cg")
    lmp.command("min_modify dmax 0.5")
    lmp.command("minimize 1.0e-5 1.0e-8 10000 100000")

lmp.command(
    "velocity all create {} {} mom yes rot yes dist gaussian".format(
        lmpsettings.tstart, np.random.randint(29847587)
    )
)

lmp.command("run {}".format(lmpsettings.runsteps))

# tidy up before closing
lmpsettings.unfix_undump(pylmp, lmp)
# lmp.command("reset_timestep 0")
lmp.command("write_restart {}".format(lmpsettings.output_lmprst))
lmp.command("clear")

# close lammps
lmp.close()
