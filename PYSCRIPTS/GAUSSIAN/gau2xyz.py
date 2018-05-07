#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import ag_unify_md as agum

"""
Convert a gaussian output-/log-file (very last entry) to a xyz-file.
"""

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    description=""
)

parser.add_argument("gau_out",
                    metavar="foo.out|foo.log",
                    help="Gaussian output-/logfile")

parser.add_argument("-out",
                    metavar="bar.xyz",
                    default="bar.xyz",
                    help="Name of xyz-file which will be written.")

args = parser.parse_args()

# read gaussian log file
gau_output  = agum.Unification()
gau_output.read_gau_log(args.gau_out)
gau_output.write_xyz(args.out)
