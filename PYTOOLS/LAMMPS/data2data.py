#!/usr/bin/env python

import argparse
import ag_unify_md as agum

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="data2data.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Convert a lammps-data-file to another with cgcmm-info.",
)

parser.add_argument(
    "-dat",
    metavar="*.lmpdat",
    action="store",
    required=True,
    help="lammps' data-file(s);\n"
    + "box-, force field-, "
    + "topology-parameters must be included!\n"
    + "http://lammps.sandia.gov/doc/2001/data_format.html",
)

parser.add_argument(
    "-out", default="DEFAULTNAME", action="store", help="Name of output-file."
)

args = parser.parse_args()

# Read Data --------------------------------------------------------------------
mydata = agum.Unification()
mydata.read_lmpdat(args.dat)
# assign atom-types if None has already been assigned
# mydata.guess_atomtypes(by_typename=True)
mydata.guess_atomtypes(by_mass=True)
# convert molecules (defined by bonds-section) to group-ids
mydata.mols_to_grps()
mydata.change_indices(incr=1, mode="increase")

# get title line
with open(args.dat, "r") as dat_in:
    mydata_title = dat_in.readline()

mydata_title = mydata_title.strip()
if "CGCMM style" in mydata_title:
    mydata_title = mydata_title.strip("CGCMM style;")

mydata.write_lmpdat(args.out + ".lmpdat", title=mydata_title, cgcmm=True)
