#!/usr/bin/env python
from __future__ import print_function
import argparse
import ag_unify_md as agum

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="placeholder.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Convert xyz-files to amber trajectory. Initially written for \
                paramfit (AmberTools) to convert xyz-files to amber-coordinate-files")

parser.add_argument("-xyz",
                    nargs="*",
                    metavar="*.xyz",
                    action="store",
                    required=True,
                    help="XYZ-File with coordinates",
                    )

parser.add_argument("-out",
                    default="DEFAULTNAME.mdcrd",
                    action="store",
                    help="Name of output-file."
                    )

args = parser.parse_args()

# parse xyz-files --------------------------------------------------------------
amber_mdcrd = agum.Unification()

# read xyz-files
for cxyz in args.xyz:
    amber_mdcrd.read_xyz(cxyz)

amber_mdcrd.write_mdcrd(args.out, range(len(amber_mdcrd.ts_coords)))
