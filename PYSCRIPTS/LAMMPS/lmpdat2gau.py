#!/usr/bin/env python

from __future__ import print_function, division
import argparse
import ag_unify_md as agum

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="lmpdat2xyz.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Convert a Lammps data file to the gaussian format."
)

parser.add_argument("lmpdat",
                    metavar="foo.lmpdat",
                    help="Lammps data file. Provides topology and force field parameters."
                    )

parser.add_argument("job_settings",
                    metavar="settings",
                    help="Settings line for gaussian."
                    )

parser.add_argument("-multiplicities",
                    nargs="*",
                    metavar="1 1 1",
                    default=[1],
                    type=int,
                    help="Multiplicit(ies) of one or several molecules/fragments"
                    )

parser.add_argument("-name_by_mass",
                    action="store_true",
                    help="Name the atom names by their types or keep the current " +
                          "naming scheme"
                    )

parser.add_argument("-sysname",
                    metavar="UNK",
                    default="UNK",
                    help="Name of the molecule/system")

parser.add_argument("-output_name",
                    metavar="foobar.gau",
                    default="foobar.gau",
                    help="Name of output-file."
                    )

args = parser.parse_args()

# File type conversion ---------------------------------------------------------

# read lammps data file
lmpdat = agum.Unification()
lmpdat.read_lmpdat(args.lmpdat)

if args.name_by_mass is True:
    lmpdat.guess_atomtypes(by_mass=True, overwrite=True)

charge = int(round(sum([i.chge for i in lmpdat.atoms])))
print(charge)
lmpdat.gaussian_charges = [charge]
lmpdat.gaussian_multiplicities = args.multiplicities
lmpdat.job_settings = args.job_settings
lmpdat.change_indices()
lmpdat.write_gau(args.output_name, -1, None, title="Converted lammps data file")
