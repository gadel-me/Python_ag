#!/usr/bin/env python


import pdb
import argparse
import ag_unify_md as agum

parser = argparse.ArgumentParser()

#parser.add_argument("-lmpdat",
#                    default=None,
#                    help="Lammps' data-file"
#                    )

parser.add_argument("pwin")
parser.add_argument("pwout")
parser.add_argument("-sf", type=float, default=1.00, help="Scaling factor for enlarging the box and its atom coordinates")
parser.add_argument("-o", default="Default")

args = parser.parse_args()
pw_sys = agum.Unification()  # coordinates

pw_sys.read_pwin(args.pwin)
pw_sys.read_pwout(args.pwout)

# enlarge coordinates by scaling factor
for idx, coord in enumerate(pw_sys.ts_coords[-1]):
    coord *= args.sf
    pw_sys.ts_coords[-1][idx] = coord

# enlarge box by 1 %
pw_sys.ts_boxes[-1].crt_a = [i * args.sf for i in pw_sys.ts_boxes[-1].crt_a]
pw_sys.ts_boxes[-1].crt_b = [i * args.sf for i in pw_sys.ts_boxes[-1].crt_b]
pw_sys.ts_boxes[-1].crt_c = [i * args.sf for i in pw_sys.ts_boxes[-1].crt_c]

# change calculation type to single point calculation
pw_sys.pw_entries["CONTROL"]["calculation"] = "scf"
print(args.o)

# write files
pw_sys.write_pwin(-1, "{}.pwin".format(args.o))
pw_sys.change_indices()
pw_sys.write_lmpdat("{}.lmpdat".format(args.o), frame_id=-1)
