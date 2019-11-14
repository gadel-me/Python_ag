#!/usr/bin/env python

import argparse
import ag_lammps
#import math
#import md_box as mdb
#import ag_unify_md as agum


# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="top_crds2lmpdat.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Merge two or more lammps data files to a new one."
)

# args for topology/forcefield and coordinates
cell  = parser.add_mutually_exclusive_group()

parser.add_argument("-lmpdats", nargs="*", required=True, default=None, metavar="foo.lmpdat", help="Lammps data files. Provides topology and force field parameters.")
parser.add_argument("-dcds", nargs="*", required=True, default=None, metavar="foo.dcd", help="Lammps dcd files. Provides coordinates only.")
parser.add_argument("-sysname", metavar="FOOBAR", default="UNK", help="Name of the molecule/system")
cell.add_argument("-cell", nargs=6, metavar=("a", "b", "c", "alpha", "beta", "gamma"), type=float, default=None, help="Box vectors a, b, c and box angles alpha, beta, gamma (in degrees) in that particular order.")
cell.add_argument("-show_lattice_cells", default=False, action="store_true")
parser.add_argument("-delete_close_atoms", default=False, action="store_true", help="BUGGY! Atoms from the system that is added to the previous one will be deleted if in close contact with existing atoms.")
parser.add_argument("-out", metavar="foobar.lmpdat", default="foobar.lmpdat", help="Name of output-file.")
args = parser.parse_args()

# Modeling ---------------------------------------------------------------------
sys_all = []

for idx, lmpdat in enumerate(args.lmpdats):

    # read dcds files if available
    try:
        curdcd = args.dcds[idx]
    except IndexError:
        curdcd = None

    cursys = ag_lammps.read_lmpdat(lmpdat, dcd=curdcd, frame_idx_start=-1, frame_idx_stop=-1)

    # resetting the pair types
    cursys.pair_types = []

    cursys.fetch_molecules_by_bonds()
    cursys.mols_to_grps()
    sys_all.append(cursys)

    # show all cell vectors of each loaded system in lattice form
    if args.show_lattice_cells is True:
        cursys.ts_boxes[-1].box_lmp2lat()
        print("System {}".format(idx))
        print(cursys.ts_boxes[-1].ltc_a)
        print(cursys.ts_boxes[-1].ltc_b)
        print(cursys.ts_boxes[-1].ltc_c)
        print(cursys.ts_boxes[-1].ltc_alpha)
        print(cursys.ts_boxes[-1].ltc_beta)
        print(cursys.ts_boxes[-1].ltc_gamma)
        print("\n")

# convert box; do it here so the linked cells are built correctly
if args.cell is not None:
    sys_all[0].ts_boxes[-1].box_lmp2lat()
    sys_all[0].ts_boxes[-1].ltc_a = args.cell[0]
    sys_all[0].ts_boxes[-1].ltc_b = args.cell[1]
    sys_all[0].ts_boxes[-1].ltc_c = args.cell[2]
    sys_all[0].ts_boxes[-1].ltc_alpha = args.cell[3]
    sys_all[0].ts_boxes[-1].ltc_beta = args.cell[4]
    sys_all[0].ts_boxes[-1].ltc_gamma = args.cell[5]
    sys_all[0].ts_boxes[-1].box_lat2lmp()

if args.show_lattice_cells is True:
    print("Exiting.")
    exit()

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

    # delete only atoms from the additional system, i.e. keep the current system
    # as is; tolerated minimal distance is 1.0
    if args.delete_close_atoms is True:
        sys_all[0].create_linked_cells(-1)
        close_atoms = sys_all[0].chk_atm_dist(min_dist=1.0)
        delete_atoms = [i for i in close_atoms if i > natoms]
        sys_all[0].delete_atoms(*delete_atoms)
        sys_all[0].refresh()
        sys_all[0].fetch_molecules_by_bonds()
        sys_all[0].mols_to_grps()
        #print("Close atoms")
        #print(close_atoms)
        #print("Delete atoms")
        #print(delete_atoms)

sys_all[0].pair_types = []
#sys_all[0].mix_pair_types(mode="ij")
#sys_all[0].fetch_molecules_by_bonds()
#sys_all[0].mols_to_grps()
sys_all[0].change_indices(incr=1, mode="increase")


# check interatomic distances fo the last time (not necessary, only for testing)
sys_all[0].create_linked_cells(-1)
close_atoms = sys_all[0].chk_atm_dist(min_dist=1.0)
close_atoms = " ".join(map(str, close_atoms))
print("***Warning: Atoms with close contacts found: \n" + close_atoms)

sys_all[0].write_lmpdat(args.out, frame_id=-1, title=args.sysname, cgcmm=True)
