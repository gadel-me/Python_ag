#!/usr/bin/env python


import argparse
import ag_unify_md as agum

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="lmpdat2xyz.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Convert a Lammps data file to the xyz format. Charges are " +
                "printed to the 5th column."
)

parser.add_argument("lmpdat",
                    metavar="foo.lmpdat",
                    help="Lammps data file. Provides topology and force field parameters."
                    )

parser.add_argument("-name_by_type",
                    action="store_true",
                    help="Name the atom names by their types or keep the current " +
                          "naming scheme"
                    )

parser.add_argument("-sysname",
                    metavar="UNK",
                    default="UNK",
                    help="Name of the molecule/system")

parser.add_argument("-out",
                    metavar="foobar.xyz",
                    default="foobar.xyz",
                    help="Name of output-file."
                    )

args = parser.parse_args()

# File type conversion ---------------------------------------------------------

# read lammps data file
lmpdat = agum.Unification()
lmpdat.read_lmpdat(args.lmpdat)

if args.name_by_type is True:
    lmpdat.guess_atomtypes(by_typename=True, overwrite=True)

lmpdat.write_xyz(args.out)
