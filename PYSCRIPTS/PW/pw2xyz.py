#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import ag_unify_md as agum
import pdb

parser = argparse.ArgumentParser()

parser.add_argument("pwscf_out",
                    metavar="*.pwscf_out",
                    help="Quantum Espresso (PW) output file (needed for coordinates and box vectors).")

parser.add_argument("-o",
                    default="foo.xyz",
                    metavar="*.xyz",
                    help="XYZ-file with coordinates and box vectors from 'pwscf_out'"
                    )

args = parser.parse_args()
pw_file_handler = agum.Unification()

# read last frame from pw output file
pw_file_handler.read_pwout(args.pwscf_out)
#pdb.set_trace()
pw_file_handler.ts_boxes[-1].box_cart2lat()
print(pw_file_handler.ts_boxes[-1].boxtype)
pw_file_handler.write_xyz(args.o, "CBZI- pwout", False, -1)
