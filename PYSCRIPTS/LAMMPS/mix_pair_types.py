#!/usr/bin/env python

import argparse
import ag_lammps

# import math
# import md_box as mdb
# import ag_unify_md as agum
import pdb

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    description="Merge two or more lammps data files to a new one.",
)

# args for topology/forcefield and coordinates
cell = parser.add_mutually_exclusive_group()

parser.add_argument(
    "lmpdats",
    nargs="*",
    metavar="foo.lmpdat",
    help="Lammps data files. Provides topology and force field parameters.",
)

parser.add_argument(
    "-out",
    metavar="foobar.lmpdat",
    default="foobar.lmpdat",
    help="Name of output-file.",
)
args = parser.parse_args()

sys_all = []

for idx, lmpdat in enumerate(args.lmpdats):
    cursys = ag_lammps.read_lmpdat(
        lmpdat, frame_idx_start=-1, frame_idx_stop=-1
    )

    # resetting the pair types
    cursys.pair_types = []

    cursys.fetch_molecules_by_bonds()
    cursys.mols_to_grps()
    sys_all.append(cursys)

# extend first molecular system with all others
for idx in range(1, len(sys_all)):

    # get number of cells from the current merged system
    if args.delete_close_atoms is True:
        natoms = len(sys_all[idx - 1].atoms)
        print(natoms)

    sys_all[0].extend_universe(sys_all[idx], u1_frame_id=-1, u2_frame_id=-1)
    sys_all[0].refresh()
    sys_all[0].fetch_molecules_by_bonds()
    sys_all[0].mols_to_grps()

#sys_all[0].pair_types = []
# ! does somehow not work properly!
# ! (old paircoeffs are written additionally to the new ones)
#pdb.set_trace()
sys_all[0].change_indices(incr=1, mode="increase")
sys_all[0].mix_pair_types(mode="ij", to_file="pair_ij.txt")
# sys_all[0].fetch_molecules_by_bonds()
# sys_all[0].mols_to_grps()

# check interatomic distances fo the last time (not necessary, only for testing)
sys_all[0].create_linked_cells(-1)
close_atoms = sys_all[0].chk_atm_dist(min_dist=1.0)
close_atoms = " ".join(map(str, close_atoms))
print("***Warning: Atoms with close contacts found: \n" + close_atoms)

#pdb.set_trace()
sys_all[0].write_lmpdat(args.out, frame_id=-1, title=args.sysname, cgcmm=True)
