#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import argparse
from lammps import lammps
from mpi4py import MPI
import os


parser = argparse.ArgumentParser()

parser.add_argument("-lmpdat",
                    default=None,
                    help="Lammps' data-file"
                    )

parser.add_argument("-lmprst",
                    default="None",
                    help="Lammps' restart file"
                    )

parser.add_argument("-set",
                    metavar="*.lmpcfg",
                    required=True,
                    help="lammps' input-file/-script with basic simulation " +
                         "settings"
                    )

parser.add_argument("-min_style",
                    default="cg",
                    help="Minimization algorithm: cg or hftn or sd or quickmin or fire. See lammps manual for more info."
                    )

parser.add_argument("-group",
                    nargs="*",
                    default="all",
                    metavar="water id < 20",
                    help="Lammps group command: ID style args. This flag may be utilized several times for several groups."
                    )

parser.add_argument("-steps",
                    metavar=500000,
                    default=500000,
                    type=int,
                    help="Number of steps after which the minimization will end (converged or not)."
                    )

parser.add_argument("-relax_box",
                    nargs=2,
                    default=None,
                    metavar="iso 1.0",
                    help="Relax box by using this given restriction (e.g. iso, tri) and pressure."
                    )

parser.add_argument("-logsteps",
                    metavar=1000,
                    default=1000,
                    type=int,
                    help="log thermodynamic-, steps- and restart-files every" +
                         "'logsteps' steps"
                    )

parser.add_argument("-out",
                    default="DEFAULTNAME",
                    help="output name for all output-files"
                    )

args = parser.parse_args()

#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator

#==============================================================================#
# Lammps
#==============================================================================#

thermargs = ["step", "temp", "press", "vol", "density",
             "cella", "cellb", "cellc", "cellalpha", "cellbeta", "cellgamma",
             "etotal", "pe", "evdwl", "ecoul", "ebond", "eangle", "edihed", "eimp",
             "enthalpy"]

# perform conjugate gradient minimization
lmp = lammps()
lmp.command("log {}.lmplog".format(args.out))
lmp.file(args.set)
lmp.command("box tilt large")  # ignore too tilted boxes

# read file
if os.path.isfile(args.lmprst) is True:
    lmp.command("read_restart {}".format(args.lmprst))
elif os.path.isfile(args.lmpdat) is True:
    lmp.command("read_data {}".format(args.lmpdat))
else:
    raise IOError("No data nor restart files given.")

# thermo output
lmp.command("thermo_style custom " + " ".join(thermargs))
lmp.command("thermo_modify lost warn flush yes")
lmp.command("thermo {}".format(args.logsteps))

# restart and dcd
lmp.command("dump dump_requench all dcd {} {}.dcd".format(args.logsteps, args.out))
lmp.command("dump_modify dump_requench unwrap yes")
lmp.command("neigh_modify every 1 delay 0 check yes")
lmp.command("restart {0} tmp-{1}.lmprst tmp-{1}.lmprst".format(args.logsteps*50, args.out))

# grouping
group = " ".join([i for i in args.group])
groupname = args.group[0].split()[0]
lmp.command("group {}".format(group))

# minimize algorithm
lmp.command("min_style {}".format(args.min_style))

if args.relax_box is not None:
    lmp.command("fix minimize_requench {} box/relax {rb[0]} {rb[1]}".format(groupname, rb=args.relax_box))

# minimize
lmp.command("minimize 1e-6 1e-8 {} 100000".format(args.steps))

# write restart file
lmp.command("write_restart {}.lmprst".format(args.out))
lmp.close()
