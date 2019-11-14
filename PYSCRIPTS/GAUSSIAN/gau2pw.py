#!/usr/bin/env python

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

parser.add_argument("pwin",
                    metavar="foo.pwscf_in",
                    help="PW input file")

parser.add_argument("-out",
                    metavar="bar.pwscf_in",
                    default="bar.pwscf_in",
                    help="Name of pw input-file which will be written.")

args = parser.parse_args()

# read gaussian log file
system  = agum.Unification()
system.read_pwin(args.pwin)
system.read_gau_log(args.gau_out, read_summary=True)
system.write_pwin(-1, args.out)
