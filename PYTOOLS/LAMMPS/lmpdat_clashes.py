#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import itertools as it
from numpy import linalg as LA
import ag_unify_md as agum

__version__ = "2017-06-02"

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(prog="refine_structure.py",
                                 formatter_class=argparse.RawTextHelpFormatter,
                                 description="Minimize energy by cooling " +
                                 "down the system to a final temperature of 5 K.")

parser.add_argument("-lmpdat",
                    metavar="*.lmpdat",
                    action="store",
                    required=True
                    )

parser.add_argument("-min_z",
                    action="store",
                    type=float,
                    required=True
                    )

parser.add_argument("-max_z",
                    action="store",
                    type=float,
                    required=True
                    )

parser.add_argument("-max_dist",
                    type=float,
                    default=1.0
                    )

parser.add_argument("-o",
                    "--outputname",
                    dest="o",
                    default="DEFAULTNAME",
                    action="store",
                    help="output file names"
                    )

args = parser.parse_args()

# Read Data --------------------------------------------------------------------
mydata = agum.Unification()
mydata.read_lmpdat(args.lmpdat)

# Find smallest distances ------------------------------------------------------
region_idx = []

# define a region
for cidx, ccoord in enumerate(mydata.ts_coords[0]):

    # z-coordinate between min- and max-z
    if args.min_z < ccoord[2] < args.max_z:
        region_idx.append(cidx)

# find atoms with small distances in that region
small_dist_idx = []

for a, b in it.combinations(region_idx, 2):
    ccoords_a = mydata.ts_coords[0][a]
    ccoords_b = mydata.ts_coords[0][b]
    # calculate distance between a and b
    vect_ab = ccoords_b-ccoords_a
    dist_ab = LA.norm(vect_ab)

    if dist_ab < args.max_dist:
        small_dist_idx.append(a)
        small_dist_idx.append(b)

small_dist_idx = set(small_dist_idx)
del region_idx

# Convert indices to ids -------------------------------------------------------
small_dist_ids = []

for cidx in small_dist_idx:
    cid = mydata.atm_idx_id[cidx]
    small_dist_ids.append(cid)

del small_dist_idx

# Get molecules ----------------------------------------------------------------
molecules = []

for cmol in mydata.molecules:
    for cid in small_dist_ids:
        if cid in cmol:
            molecules.extend(cmol)
            break

print(molecules)
#molecules = [i-1 for i in molecules]
#molecules = [str(i) for i in molecules]
#molecules = " ".join(molecules)
#print(molecules)
