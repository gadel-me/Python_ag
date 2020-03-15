#!/usr/bin/env python


import math
import numpy as np
import argparse
import md_box as mdb
import ag_unify_md as agum
import ag_vectalg as agv

"""
Cut a shape by its box vectors from a given lammps file. The shape or the negative
of it may be returned. Cell from lammps data file must be unwrapped.qweqwee
"""

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("lmpdat", metavar="*.lmpdat", help="lammps' data-file")

parser.add_argument("-dcd", default=None, help="DCD file for additional coordinates")

parser.add_argument("-frame", type=int, default=0, help="Frame to cut coordinates from")

parser.add_argument(
    "-shift_main",
    nargs=3,
    type=float,
    default=(0.0, 0.0, 0.0),
    metavar=("sx, sy, sz"),
    help="Shift lmpdat by vector sx sy sz. TBD",
)

parser.add_argument(
    "-shift_main_by_cog",
    default=False,
    action="store_true",
    help="Shift lmpdats coordinates so that its center of mass corresponds P(0/0/0) (coordinates from dcd will be shifted as well).",
)

parser.add_argument(
    "-shift",
    nargs=3,
    type=float,
    default=(0.0, 0.0, 0.0),
    metavar=("px, py, pz"),
    help="Cartesian coordinates of the positional-vector which shifts the box.",
)

parser.add_argument(
    "-box",
    required=True,
    type=float,
    nargs="*",
    help="Box vectors in lammps, cartesian or lattice.",
)

parser.add_argument(
    "-boxtype",
    required=True,
    help="Box types: 'lammps' | 'cartesian' | 'lattice'\n"
    + "Lattice: a, b, c, alpha, beta, gamma in angstrom/degrees\n"
    + "Lammps: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz",
)

parser.add_argument(
    "-lmpbox",
    default=False,
    metavar="*.lmpdat",
    help="lammps' data-file with box coordinates to cut.",
)

parser.add_argument(
    "-dcdbox",
    default=False,
    metavar="*.lmpdat",
    help="lammps' dcd-file with box coordinates to cut.",
)

parser.add_argument(
    "-dcdbox_frame",
    type=int,
    default=-1,
    metavar="*.lmpdat",
    help="frame of the dcd file to use box information from.",
)

parser.add_argument(
    "-inverse",
    default=False,
    action="store_true",
    help="Selection: box (normal), negative of box (inverted)",
)

args = parser.parse_args()

sys = agum.Unification()
print("Reading...")
sys.read_lmpdat(args.lmpdat)

if args.dcd is not None:
    sys.import_dcd(args.dcd)
    sys.read_frames(frame=args.frame, to_frame=args.frame + 1)

# sys_box = copy.deepcopy(sys.ts_boxes[-1])
# sys_box.box_lmp2cart()

if args.lmpbox is False or args.dcdbox is False:
    if args.boxtype == "lammps":
        box = mdb.Box(
            boxtype=args.boxtype,
            lmp_xlo=args.box[0],
            lmp_xhi=args.box[1],
            lmp_ylo=args.box[2],
            lmp_yhi=args.box[3],
            lmp_zlo=args.box[4],
            lmp_zhi=args.box[5],
            lmp_xy=args.box[6],
            lmp_xz=args.box[7],
            lmp_yz=args.box[8],
        )
        box.box_lmp2cart()  # change box vectors to lattice
    elif args.boxtype == "lattice":
        box = mdb.Box(
            boxtype=args.boxtype,
            ltc_a=args.box[0],
            ltc_b=args.box[1],
            ltc_c=args.box[2],
            ltc_alpha=math.radians(args.box[3]),
            ltc_beta=math.radians(args.box[4]),
            ltc_gamma=math.radians(args.box[5]),
        )
        box.box_lat2cart()
    else:
        box = mdb.Box(
            boxtype=args.boxtype,
            crt_a=[args.box[0], args.box[1], args.box[2]],
            crt_b=[args.box[3], args.box[4], args.box[5]],
            crt_c=[args.box[6], args.box[7], args.box[8]],
        )
else:
    sys_box = agum.Unification()

    if args.lmpbox:
        sys_box.read_lmpdat(args.lmpbox)

# define box faces
face_1 = agv.get_plane(box.crt_a, box.crt_b, vt_p=args.shift)  # xy
face_2 = agv.get_plane(box.crt_c, box.crt_a, vt_p=args.shift)  # xz
face_3 = agv.get_plane(box.crt_b, box.crt_c, vt_p=args.shift)  # yz
# multiply with -1 because normal vector must show in the opposite direction
face_4 = [
    -1 * i
    for i in agv.get_plane(
        box.crt_a, box.crt_b, vt_p=agv.add_vts(box.crt_c, args.shift)
    )
]
face_5 = [
    -1 * i
    for i in agv.get_plane(
        box.crt_c, box.crt_a, vt_p=agv.add_vts(box.crt_b, args.shift)
    )
]
face_6 = [
    -1 * i
    for i in agv.get_plane(
        box.crt_b, box.crt_c, vt_p=agv.add_vts(box.crt_a, args.shift)
    )
]

pdb.set_trace()

if args.shift_main_by_cog is True:
    sys.transpose_by_cog(args.frame, np.array([0, 0, 0]), copy=False)

print("Cutting shape...")
sys.cut_shape(-1, args.inverse, face_1, face_2, face_3, face_4, face_5, face_6)

# change box if cut is inverse
if args.inverse is True:
    box.box_cart2lmp()
    sys.ts_boxes[-1] = box

sys.refresh()
sys.change_indices(incr=1, mode="increase")
print("Writing lammps data file...")
sys.write_lmpdat(
    "{}_out.lmpdat".format(args.lmpdat.rstrip(".lmpdat")),
    frame_id=-1,
    title="cut box",
    cgcmm=True,
)
