#!/usr/bin/env python
from __future__ import print_function, division
import argparse
#import os
import ag_gaussian as aggau
from natsort import natsorted

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    description="Extract energy from gaussian log-/output-file.")

parser.add_argument("gau_logfiles",
                    nargs="*",
                    metavar="*.log/*.out",
                    help="Gaussian log-/output-file. Works also with '*' if all files in directory should be considered (only files ending with '.log' will be considered)."
                    )

parser.add_argument("-out",
                    metavar="*.txt",
                    help="File in which the energy values are written.")

args = parser.parse_args()

# process gaussian output ------------------------------------------------------
gau = aggau.GauStuff()
# only use files ending with ".log", get them in the right order
logfiles = natsorted([i for i in args.gau_logfiles if i.endswith(".log")])

# extract energy/ies from each file
for cur_log in logfiles:
    gau.read_log(cur_log)

# write simple output (e.g. for use with gnuplot)
with open(args.out, "w") as output:
    for energy in gau.gau_energies:
        output.write("{}\n".format(energy))
