#!/usr/bin/env python

import argparse
#import math
#import md_box as mdb
import ag_unify_md as agum

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="top_crds2lmpdat.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Replicate a lammps data file in a, b and c direction to form \
    a super cell."
)

parser.add_argument("lmpdat",
                    metavar="foo.lmpdat",
                    help="Lammps data file. Provides topology, force field parameters and cell vectors."
                    )

parser.add_argument("-dcd",
                    default=None,
                    metavar="foo.dcd",
                    help="CHARMM dcd file."
                    )

parser.add_argument("-frame",
                    type=int,
                    default=0,
                    metavar="0",
                    help="Number of frame to read. Caveat: '-2' if last frame is desired."
                    )

parser.add_argument("-a",
                    nargs=2,
                    type=int,
                    default=(0, 0),
                    metavar=12,
                    help="replicate cell in cell vector a direction this many times.")

parser.add_argument("-b",
                    nargs=2,
                    type=int,
                    default=(0, 0),
                    metavar=5,
                    help="replicate cell in cell vector b direction this many times.")

parser.add_argument("-c",
                    nargs=2,
                    type=int,
                    default=(0, 0),
                    metavar=8,
                    help="replicate cell in cell vector c direction this many times.")

parser.add_argument("-sysname",
                    metavar="FOOBAR",
                    default="UNK",
                    help="Name of the molecule/system")

parser.add_argument("-mix_pair_types",
                    action="store_true",
                    help="Mix pair types if Pair Coeff given in data file")

parser.add_argument("-out",
                    metavar="foobar.lmpdat",
                    default="foobar.lmpdat",
                    help="Name of output-file."
                    )

#parser.add_argument("-cgcmm",
#                    action="store_true",
#                    help="Write file with cgcmm style info.")

parser.add_argument("-guess_atoms",
                    action="store_true",
                    help="Guess atoms by mass.")

args = parser.parse_args()

# Modeling ---------------------------------------------------------------------
sys  = agum.Unification()
sys.read_lmpdat(args.lmpdat, energy_unit="eV", angle_unit="deg")

if args.guess_atoms is True:
    sys.guess_atomtypes(by_mass=True)

if args.dcd is not None:
    sys.import_dcd(args.dcd)
    sys.read_frames(frame=args.frame, to_frame=args.frame + 1)
    sys.close_dcd()

sys.replicate_cell(n_start=args.a[0], n_stop=args.a[1], direction="a", frame_id=-1, adjust_box=True)
sys.replicate_cell(n_start=args.b[0], n_stop=args.b[1], direction="b", frame_id=-1, adjust_box=True)
sys.replicate_cell(n_start=args.c[0], n_stop=args.c[1], direction="c", frame_id=-1, adjust_box=True)
sys.fetch_molecules_by_bonds()
sys.mols_to_grps()

if args.mix_pair_types is True:
    sys.mix_pair_types(mode="ij")

sys.change_indices(incr=1, mode="increase")
sys.write_lmpdat(args.out, frame_id=-1, title=args.sysname, cgcmm=True)
