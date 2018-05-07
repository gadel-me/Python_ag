#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import ag_unify_md as agum

parser = argparse.ArgumentParser()

parser.add_argument("pwscf_in",
                    metavar="*.pwscf_in",
                    help="Quantum Espresso (PW) input file (needed for settings)."
                    )

parser.add_argument("pwscf_out",
                    metavar="*.pwscf_in",
                    help="Quantum Espresso (PW) output file (needed for coordinates and box vectors)."
                    )

parser.add_argument("-frame_id",
                    type=int,
                    default=-1,
                    metavar="-1",
                    help="ID of Frame from pwscf_out to use coordinates and cell from."
                    )

parser.add_argument("-o",
                    default="foo.pwscf_in",
                    metavar="*.pwscf_in",
                    help="Input file with coordinates and box vectors from 'pwscf_out'"
                    )

args = parser.parse_args()
pw_file_handler = agum.Unification()
pw_file_handler.read_pwin(args.pwscf_in)
pw_file_handler.read_pwout(args.pwscf_out)
pw_file_handler.write_pwin(-1, args.o)
