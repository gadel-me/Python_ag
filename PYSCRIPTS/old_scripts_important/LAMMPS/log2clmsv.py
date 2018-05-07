#!/usr/bin/env python
from __future__ import print_function, division
#import readline  # necessary for raw_input and using arrow keys
import argparse
import ag_unify_log as agul

__version__ = "2017-05-26"

# Extract the thermo output from a lammps-log file and write the values in
# a nicer format (clmsv) which can be plotted more easily by gnuplot,
# xmgrace, etc.

parser = argparse.ArgumentParser(prog="log2clmsv.py",
                                 formatter_class=argparse.RawTextHelpFormatter,
                                 description="Read one or more log.lammps and " +
                                             "convert them to clmsv files. " +
                                             "Entries are split to enumerated files.")

parser.add_argument("-log",
                    metavar="log.lammps",
                    action="store",
                    required=True,
                    nargs="*",
                    help="Column separated files containing values " +
                         "for several keywords"
                    )

parser.add_argument("-out",
                    action="store",
                    required=True,
                    help="Output name for clmsv-files."
                    )

parser.add_argument("-s",
                    "--split",
                    action="store_true",
                    help="If log.lammps consists of several entries, " +
                         "split output into multiple clmsv-files.")

args = parser.parse_args()

a = agul.LogUnification()
a.read_lmplog(*args.log)
a.write_clmsv(args.out, split=args.split)
