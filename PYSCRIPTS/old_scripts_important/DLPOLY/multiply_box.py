#!/home/angad/Programs/anaconda/bin/python -i

from __future__ import print_function, division
#import os
#import tqdm
#import numpy as np
import sys
sys.path.append("/home/gadelmeier/Python/myPYMODULES/DLPOPLY_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/OTHER_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/MATH_MODULES/")
import DL_HISTORY_class8_new_approaches_8_4 as dlhc
#import my_lin_alg_compendium as mlc
import argparse

program_version = "Multiply CONFIG-box, build 0.01"

# Argument parsing stuff *******************************************************
parser = argparse.ArgumentParser(prog="multiply_box.py",
                                 description="Multiply a given box by" +
                                 "its cell vectors")

parser.add_argument("-cin",
                    "--dlpolyconfig",
                    dest="dlpolycfg_in",
                    metavar="xy.dlpolycfg",
                    default=None,
                    action="store",
                    help="DL_POLY-CONFIG- or -REVCON-File. " +
                    "Box vectors must be present in File!",
                    )

parser.add_argument("-cout",
                    "--dlpolyconfig_out",
                    dest="dlpolycfg_out",
                    metavar="out.dlpolycfg",
                    default=None,
                    action="store",
                    help="Filename of modified xy.dlpolycfg." +
                    "If no 'out.dlpolycfg' is defined, no file will be written",
                    )

parser.add_argument("-fld",
                    "--dlpolyfld",
                    metavar="xy.dlpolyfld",
                    action="store",
                    required=True,
                    help="DL_POLY-FIELD-File, belonging " +
                         "to xy.dlpolycfg/xy.dlpolyhistory",
                    )

parser.add_argument("-a",
                    "--box_vector_a",
                    dest="num_a_additions",
                    metavar="Box Vector 'a'",
                    default=None,
                    type=int,
                    action="store",
                    help="Times how often vector a will be added to existing" +
                    "vector a",
                    )

parser.add_argument("-b",
                    "--box_vector_b",
                    dest="num_b_additions",
                    metavar="Box Vector 'b'",
                    default=None,
                    type=int,
                    action="store",
                    help="Times how often vector b will be added to existing" +
                    "vector b",
                    )

parser.add_argument("-c",
                    "--box_vector_c",
                    dest="num_c_additions",
                    metavar="Box Vector 'c'",
                    default=None,
                    type=int,
                    action="store",
                    help="Times how often vector c will be added to existing" +
                    "vector c",
                    )

args = parser.parse_args()

# PROCESS FIELD ****************************************************************
finfo = dlhc.FieldFile(args.dlpolyfld)
finfo.read_field()

# PROCESS CONFIG ***************************************************************
conf = dlhc.ConfigFile(finfo, args.dlpolycfg_in)
conf.read_config()
# Get number of total atoms in CONFIG
conf.config.refresh_natms()


if args.num_a_additions:
    conf.config.multiply_config(box_vector="a", n=args.num_a_additions)
    conf.config.refresh_natms()
if args.num_b_additions:
    conf.config.multiply_config(box_vector="b", n=args.num_b_additions)
    conf.config.refresh_natms()
if args.num_c_additions:
    conf.config.multiply_config(box_vector="c", n=args.num_c_additions)
    conf.config.refresh_natms()

conf.config.write_config(config_out=args.dlpolycfg_out,
                         file_header=conf.config_header)
