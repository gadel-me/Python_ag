#!/usr/bin/env python

import argparse
import itertools as it
from numpy import linalg as LA
import pdb
import ag_unify_md as agum

__version__ = "2017-06-02"

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="", formatter_class=argparse.RawTextHelpFormatter, description=""
)

parser.add_argument("-lmpdat", metavar="*.lmpdat", action="store", required=True)

parser.add_argument(
    "-dcd", required=True, metavar="*.dcd", action="store", help="Lammps' DCD-file."
)

parser.add_argument(
    "-f",
    "--frame",
    default=-1,
    type=int,
    help="Index (!) of frame to convert (negative indices allowed).",
)

parser.add_argument("-min_dist", type=float, default=1.0)


args = parser.parse_args()

# Read Data --------------------------------------------------------------------
mydata = agum.Unification()

frame = args.frame
mydata = agum.Unification()
mydata.read_lmpdat(args.lmpdat)
mydata.import_dcd(args.dcd)
mydata.read_frames(frame=None, to_frame=-1, frame_by="index")
mydata.create_linked_cells(frame_id=args.frame, rcut_a=2, rcut_b=2, rcut_c=2)
close_contacts = mydata.chk_atm_dist(
    frame_id=args.frame, min_dist=args.min_dist, exclude_same_molecule=False
)
# pdb.set_trace()
print(close_contacts)
