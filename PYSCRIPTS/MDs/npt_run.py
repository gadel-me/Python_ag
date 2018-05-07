#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import argparse
import os
import re
from lammps import lammps
from mpi4py import MPI

"""
Caveat: Ice cube prevention and velocity creation is used on all atoms.
"""

parser = argparse.ArgumentParser()

parser.add_argument("-lmpdat",
                    default=None,
                    help="Lammps' data-file"
                    )

parser.add_argument("-lmprst",
                    default=None,
                    help="Lammps' restart file"
                    )

parser.add_argument("-gpu",
                    nargs=2,
                    metavar=("neigh", "no"),
                    default=None,
                    help="Utilize gpu-package and build neighbor (neigh) lists on the gpu (yes) or on the cpu (no).",
                    )

parser.add_argument("-set",
                    metavar="*.lmpcfg",
                    required=True,
                    help="lammps' input-file/-script with basic simulation " +
                         "settings"
                    )

parser.add_argument("-group",
                    nargs="*",
                    default="all",
                    metavar="water id < 20",
                    help="Lammps group command: ID style args. This flag may be utilized several times for several groups."
                    )

parser.add_argument("-ensemble",
                    nargs="*",
                    default=None,
                    metavar="FIXNPTISO npt temp 150 290 0.1 iso 1.0 1.0 1",
                    help="Hoover npt settings for lammps. Currently only hoover " +
                         "ensembles can be utilized. If not chosen, no md will be run."
                    )

parser.add_argument("-minimize",
                    action="store_true",
                    help="Do a minimization before the md run."
                    )

parser.add_argument("-min_style",
                    default="cg",
                    help="Minimization algorithm: cg or hftn or sd or quickmin or fire. See lammps manual for more info."
                    )

parser.add_argument("-relax_box",
                    nargs=2,
                    default=None,
                    metavar="iso 1.0",
                    help="Relax box by using this given restriction (e.g. iso, tri) and pressure."
                    )

parser.add_argument("-steps",
                    metavar=500000,
                    default=500000,
                    type=int,
                    help="Number of steps after which the minimization will end (converged or not)."
                    )

parser.add_argument("-logsteps",
                    metavar=1000,
                    default=1000,
                    type=int,
                    help="log thermodynamic-, steps- and restart-files every" +
                         "'logsteps' steps"
                    )

parser.add_argument("-timeout",
                    metavar="23:50:00",
                    default="23:50:00",
                    help="allowed duration of simulation, for resubmitting purposes;  should be < 24h")

parser.add_argument("-out",
                    default="DEFAULTNAME",
                    help="output name for all output-files"
                    )

parser.add_argument("-debug",
                    action="store_true",
                    help="Write all system info to a file called 'info.txt'"
                    )

args = parser.parse_args()

if args.ensemble is not None:
    start_temp = re.findall(r'temp \d*', args.ensemble[0])[0].split()[-1]

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
#lmp.command("timer timeout {} every {}".format(args.timeout, args.logsteps))

# basic lammps settings
if args.gpu is not None:
    lmp.command("package gpu 1 {neigh[0]} {neigh[1]}".format(neigh=args.gpu))
    lmp.command("suffix gpu")

lmp.file(args.set)
lmp.command("box tilt large")  # ignore too tilted boxes

# read file
if args.lmprst is not None and os.path.isfile(args.lmprst):
    lmp.command("read_restart {}".format(args.lmprst))
elif os.path.isfile(args.lmpdat) is True:
    lmp.command("read_data {}".format(args.lmpdat))

    if not args.minimize:
        lmp.command(("velocity all create {} 483806 rot yes dist gaussian").format(start_temp))

else:
    raise IOError("No data nor restart files given.")

# thermo stuff
lmp.command("thermo_style custom " + " ".join(thermargs))
lmp.command("thermo_modify lost warn flush yes")
lmp.command("thermo {}".format(args.logsteps))

# trajectory stuff
lmp.command("dump dump_requench all dcd {} {}.dcd".format(args.logsteps, args.out))
lmp.command("dump_modify dump_requench unwrap yes")

# restart stuff
#lmp.command("neigh_modify every 1 delay 0 check yes")
lmp.command("restart {0} tmp-{1}.lmprst tmp-{1}.lmprst".format(args.logsteps*50, args.out))

# ice cube prevention
lmp.command("fix ic_prevention all momentum {} linear 1 1 1 angular rescale".format(100))

# grouping
if args.group == "all":
    groupname = args.group
else:
    group = " ".join([i for i in args.group])
    groupname = args.group[0].split()[0]
    lmp.command("group {}".format(group))

# minimization stuff
if args.minimize is True:
    # choose minimize algorithm, e.g. steepest descent or conjugate gradient
    lmp.command("min_style {}".format(args.min_style))

    if args.relax_box is not None:
        lmp.command("fix minimize_requench {} box/relax {rb[0]} {rb[1]}".format(groupname, rb=args.relax_box))

    # minimize
    lmp.command("minimize 1e-6 1e-8 {} 100000".format(args.steps))

# md stuff
if args.ensemble is not None:
    args.ensemble = args.ensemble[0]
    lmp.command("fix ENSEMBLE {0} {1}".format(groupname, args.ensemble))
    lmp.command("run {}".format(args.steps))

lmp.command("write_restart {}.lmprst".format(args.out))

if args.debug is True:
    lmp.command("info all out append {}_info.txt".format(args.out))

lmp.close()
