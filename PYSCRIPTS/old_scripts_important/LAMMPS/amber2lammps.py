#!/usr/bin/env python
from __future__ import print_function
import argparse
import math
import md_box as mdb
import ag_unify_md as agum

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="data2data.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Convert an Amber Topology and a coordinate File (xyz, amber-inpcrd) to " +
                "a LAMMPS Data File")

parser.add_argument("-prmtop",
                    required=True,
                    metavar="*.prmtop",
                    action="store",
                    help="Amber Topology File",
                    )

parser.add_argument("-xyz",
                    metavar="*.xyz",
                    action="store",
                    help="XYZ-File with coordinates",
                    )

parser.add_argument("-inpcrd",
                    metavar="*.inpcrd",
                    action="store",
                    help="Amber-Coordinates file.")

parser.add_argument("-numunits",
                    type=int,
                    action='store',
                    default=1,
                    help="Times to replicate the topology to match the number \
                         of atoms/molecules in the XYZ-File")

parser.add_argument("-box",
                    nargs="*",
                    metavar="a, b, c, alpha, beta, gamma",
                    type=float,
                    default=None,
                    help="Box vectors a, b, c and box angles alpha, beta gamma\
                          in that particular order! Box angles in degrees!")

parser.add_argument("-sysname",
                    default="UNK",
                    type=str,
                    help="Name of the molecule/system")

parser.add_argument("-build_supercell",
                    nargs="*",
                    type=int,
                    default=[1, 1, 1],
                    help="Build a super cell which has n-times the size of all \
                    three directional vectors a b c")

parser.add_argument("-out",
                    default="DEFAULTNAME.lmpdat",
                    action="store",
                    help="Name of output-file."
                    )

args = parser.parse_args()

# Read Prmtop ------------------------------------------------------------------
amber_topology = agum.Unification()
amber_topology.read_prmtop(args.prmtop)
#amber_topology.mix_pair_types(mode="ij")
amber_topology.ui_convert_units(energy_unit_out='eV',
                                ang_unit_out="deg",
                                cvff_style=False)
amber_topology.mix_pair_types(mode="ij")
# replicate topology for each unit to come (e.g. prmtop has 1 molecule -> replicate n-molecules times)
amber_topology.add_topology_replicate(args.numunits)

# get coordinates
if args.xyz:
    amber_topology.read_xyz(args.xyz, overwrite_data=False)
else:
    amber_topology.read_inpcrd(args.inpcrd)

# Box definition ---------------------------------------------------------------
if args.box is None:
    amber_topology.def_boxes_by_coords()
else:
    a, b, c = args.box[0:3]
    alpha, beta, gamma = [math.radians(i) for i in args.box[3:]]
    amber_topology.ts_boxes.append(mdb.Box(boxtype="lattice",
                                   ltc_a=a, ltc_b=b, ltc_c=c,
                                   ltc_alpha=alpha, ltc_beta=beta,
                                   ltc_gamma=gamma))

# build the super cell
if args.build_supercell != [1, 1, 1]:
    for n, aim in zip(args.build_supercell, ("a", "b", "c")):
        amber_topology.replicate_cell(n, aim)

amber_topology.change_indices(incr=1, mode="increase")

# write title containing box info (if there)
if args.box is None:
    amber_topology.write_lmpdat(args.out,
                                frame_id=0,
                                title="{}".format(args.sysname),
                                cgcmm=True)
else:
    amber_topology.write_lmpdat(args.out,
                                frame_id=0,
                                title="{}, Unit-Cell: a={box[0]}, b={box[1]}, c={box[2]}".format(args.sysname, box=args.box[0:3]),
                                cgcmm=True)
