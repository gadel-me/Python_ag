#!/usr/bin/env python

import argparse
import ag_unify_md as agum

"""
Convert a gaussian output-/log-file (very last entry) to a xyz-file.
"""

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter, description=""
)

parser.add_argument(
    "gau_out", metavar="foo.out|foo.log", help="Gaussian output-/logfile"
)

parser.add_argument(
    "-out",
    metavar="bar.mdcrd",
    default="bar.mdcrd",
    help="Name of mdcrd-file which will be written.",
)

args = parser.parse_args()

# read gaussian log file
gau_output = agum.Unification()
gau_output.read_gau_log(args.gau_out)
gau_output.write_mdcrd(args.out, *list(range(len(gau_output.ts_coords))))
