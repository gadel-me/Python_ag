#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import ag_unify_md as agum

parser = argparse.ArgumentParser()

parser.add_argument("pwscf_in",
                    metavar="*.pwscf_in",
                    help="Quantum Espresso (PW) input file (needed for settings)."
                    )

parser.add_argument("-pwscf_out",
                    default=None,
                    metavar="*.pwscf_out",
                    help="Quantum Espresso (PW) output file (needed for coordinates and box vectors)."
                    )

parser.add_argument("-ecutwfc",
                    type=int,
                    default=None,
                    metavar="47",
                    help="kinetic energy cutoff (Ry) for wavefunctions"
                    )

parser.add_argument("-ecutrho",
                    type=int,
                    metavar="323",
                    default=None,
                    help="Kinetic energy cutoff (Ry) for charge density and potential"
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

# read last frame from pw output file
if args.pwscf_out is not None:
    pw_file_handler.read_pwout(args.pwscf_out)

# define ecutwfc by user input
if args.ecutwfc is not None:
    pw_file_handler.pw_entries["SYSTEM"]["ecutwfc"] = args.ecutwfc

# define ecutrho by user input
if args.ecutrho is not None:
    pw_file_handler.pw_entries["SYSTEM"]["ecutrho"] = args.ecutrho

pw_file_handler.write_pwin(args.frame_id, args.o)