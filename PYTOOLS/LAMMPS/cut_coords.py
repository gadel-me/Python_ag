#!/usr/bin/env python


import math
import argparse
import md_box as mdb
import ag_unify_md as agum
import ag_vectalg as agv
import copy
import numpy as np

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("mainfile",
                    metavar="*.lmpdat",
                    help="lammps' data-file"
                    )

parser.add_argument("coordinates_file",
                    metavar="*.xyz",
                    help="File with coordinates which are used to cut from mainfile."
                    )

parser.add_argument("-out",
                    default="foo.lmpdat",
                    metavar="foo.lmpdat",
                    help="Output-File."
                    )

args = parser.parse_args()

mainsys = agum.Unification()
mainsys.read_lmpdat(args.mainfile)
#mainsys.ts_coords[-1] = [round(i, 3) for j in mainsys.ts_coords[-1] for i in j]
cutsys = agum.Unification()
#cutsys.read_pdb(args.coordinates_file)
cutsys.read_xyz(args.coordinates_file)
delete_atoms = []

mainsys.ts_coords[-1] = [list([round(j, 3) for j in i]) for i in mainsys.ts_coords[-1]]
cutsys.ts_coords[-1] = [list([round(j, 3) for j in i]) for i in cutsys.ts_coords[-1]]

for atm_idx, mainsys_coords in enumerate(mainsys.ts_coords[-1]):
    if mainsys_coords not in cutsys.ts_coords[-1]:
        delete_atoms.append(atm_idx)

mainsys.delete_atoms(*delete_atoms)
mainsys.refresh()

delete_atoms2 = []

for mol in mainsys.molecules:
    if len(mol) == 1:
        idx = list(mol)
        delete_atoms2.append(idx[0])

mainsys.delete_atoms(*delete_atoms2)
mainsys.refresh()
mainsys.change_indices(incr=1, mode="increase")
mainsys.write_lmpdat(args.out, frame_id=-1, title="cut crystal", cgcmm=True)
