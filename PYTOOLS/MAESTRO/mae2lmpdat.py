#!/usr/bin/env python

import os
import argparse
import md_box as mdb
import ag_unify_md as agum

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
group = parser.add_mutually_exclusive_group()
parser.add_argument("mae",
                    metavar="*.mae",
                    help="Maestro-file.",
                    )

group.add_argument("-prmtop",
                   metavar="*.prmtop",
                   help="Amber topology file (for force field and topology).",
                   )

group.add_argument("-lmpdat",
                   metavar="*.lmpdat",
                   help="Lammps data file (for force field and topology).",
                   )

args = parser.parse_args()

# Read mae-file ----------------------------------------------------------------
mae = agum.Unification()

# force field from prmtop or lammps data
if args.prmtop:
    mae.read_prmtop(args.prmtop)
    mae.ui_convert_units(energy_unit_out='eV',
                         ang_unit_out="deg",
                         cvff_style=False)
    mae.mix_pair_types(mode="ij")
elif args.lmpdat:
    mae.read_lmpdat(args.lmpdat,
                    energy_unit="eV",
                    angle_unit="deg",
                    overwrite_data=True)
else:
    raise IOError("Please choose a prmtop or lmpdat file!")

# delete given coordinates and boxes
mae.ts_coords = []
mae.ts_boxes  = []

# coordinates from jaguar
mae.read_mae(args.mae)
# increment indices by 1 (internally everything starts with 0)
mae.change_indices(incr=1, mode="increase")
# trim to base name without mae ending
mae_basename = os.path.basename(args.mae).rstrip(".01.mae")

nframes = len(mae.ts_coords)  # write one file for each conformation
for cf in range(nframes):
    # create new boxes for each conformation
    cur_box = mdb.Box(boxtype="lammps",
                      lmp_xlo=-20, lmp_xhi=20,
                      lmp_ylo=-20, lmp_yhi=20,
                      lmp_zlo=-20, lmp_zhi=20)
    mae.ts_boxes.append(cur_box)
    out = "{}_{}.lmpdat".format(mae_basename, cf)
    mae.write_lmpdat(out, frame_id=cf, title=mae_basename, cgcmm=True)
