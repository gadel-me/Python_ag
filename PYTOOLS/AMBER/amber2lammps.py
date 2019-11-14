#!/usr/bin/env python

import argparse
import math
import md_box as mdb
import ag_unify_md as agum

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="amber2lammps.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Convert prmtop- and inpcrd(-or a xyz)-files from the AMBER Program to lmpdat (LAMMPS data file).")

parser.add_argument("prmtop",
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

parser.add_argument("-mdcrd",
                    metavar="*.mdcrd",
                    help="Amber-Trajectory file.")

parser.add_argument("-box",
                    nargs="*",
                    metavar="a, b, c, alpha, beta, gamma",
                    type=float,
                    default=None,
                    help="Box vectors a, b, c and box angles (in degrees) alpha, beta gamma\
                          in that particular order")

parser.add_argument("-sysname",
                    default="UNK",
                    type=str,
                    help="Name of the molecule/system")

parser.add_argument("-out",
                    default="DEFAULTNAME.lmpdat",
                    action="store",
                    help="Name of output-file."
                    )

args = parser.parse_args()

# Topology and Force Field -----------------------------------------------------
amber_topology = agum.Unification()
amber_topology.read_prmtop(args.prmtop)
# convert units from kcal to eV
amber_topology.ui_convert_units(energy_unit_out='eV',
                                ang_unit_out="deg",
                                cvff_style=False)
amber_topology.mix_pair_types(mode="ij")

# Coordinates ------------------------------------------------------------------
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

# change indices to start with 1 (everything is newly ordered)
amber_topology.change_indices(incr=1, mode="increase")

# write title containing box info (if there)
amber_topology.write_lmpdat(args.out,
                            frame_id=0,
                            title="{}".format(args.sysname),
                            cgcmm=True)
