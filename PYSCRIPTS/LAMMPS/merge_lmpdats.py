#!/usr/bin/env python
from __future__ import print_function, division
import argparse
#import math
#import md_box as mdb
import ag_unify_md as agum

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="top_crds2lmpdat.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Merge two or more lammps data files to a new one."
)

# args for topology/forcefield and coordinates
cell  = parser.add_mutually_exclusive_group()

parser.add_argument("-lmpdats",
                    nargs="*",
                    required=True,
                    default=None,
                    metavar="foo.lmpdat",
                    help="Lammps data file. Provides topology and force field parameters."
                    )

parser.add_argument("-sysname",
                    metavar="FOOBAR",
                    default="UNK",
                    help="Name of the molecule/system")

cell.add_argument("-cell",
                  nargs=6,
                  metavar=("a", "b", "c", "alpha", "beta", "gamma"),
                  type=float,
                  default=None,
                  help="Box vectors a, b, c and box angles alpha, beta, gamma" +
                       "(in degrees) in that particular order."
                  )

#cell.add_argument("-cbg",
#                  action="store_true",
#                  help="Guess cell by atomic coordinates."
#                  )
#
#cell.add_argument("-cbc",
#                  metavar="foo.lmpdat",
#                  default=None,
#                  help="Name of data file from which the coordinates should be taken."
#                  )

parser.add_argument("-out",
                    metavar="foobar.lmpdat",
                    default="foobar.lmpdat",
                    help="Name of output-file."
                    )

args = parser.parse_args()

# Modeling ---------------------------------------------------------------------
sys_all = []

for lmpdat in args.lmpdats:
    sys  = agum.Unification()
    sys.read_lmpdat(lmpdat, energy_unit="eV", angle_unit="deg")
    sys_all.append(sys)

for idx in xrange(1, len(sys_all)):
    sys_all[0].extend_universe(sys_all[idx])

sys_all[0].fetch_molecules_by_bonds()
sys_all[0].mols_to_grps()
sys_all[0].change_indices(incr=1, mode="increase")
sys_all[0].write_lmpdat(args.out, title=args.sysname, cgcmm=True)
