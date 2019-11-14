#!/usr/bin/env python

import argparse
import ag_unify_md as agum

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    description=""
)

parser.add_argument("gau_out",
                    metavar="foo.out|foo.log",
                    help="Gaussian output-/logfile"
                    )

args = parser.parse_args()

# read gaussian log file
gau_output  = agum.Unification()
gau_output.read_gau_log(args.gau_out)

# vector between cbz1-amide-hydrogen and cbz2-amide-oxygen
