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
import my_linalg3 as ml
import argparse
#import copy


def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out


parser = argparse.ArgumentParser(prog="add_box2box.py",
                                 description="Cut a new cell from a given " +
                                             "cell in DLPOLY-CONFIG-Format.")

parser.add_argument("-cin",
                    "--config_in",
                    dest="cin",
                    metavar="xy.dlpolycfg",
                    default=None,
                    action="store",
                    help="DL_POLY-CONFIG- or -REVCON-File to add cell to. " +
                    "Cell vectors must be present in File and File must be " +
                    "unwrapped!",
                    )

parser.add_argument("-cout",
                    "--config_out",
                    dest="cout",
                    metavar="out.dlpolycfg",
                    default=None,
                    action="store",
                    help="Filename of modified xy.dlpolycfg." +
                    "If no 'out.dlpolycfg' is defined, no file will be written",
                    )

parser.add_argument("-fin",
                    "--dlpolyfld",
                    dest="fin",
                    metavar="xy.dlpolyfld",
                    action="store",
                    required=True,
                    help="DL_POLY-FIELD-File, belonging " +
                         "to xy.dlpolycfg/xy.dlpolyhistory",
                    )

parser.add_argument("--vectors",
                    nargs="*",
                    metavar="direction",
                    type=float,
                    help="Direction to which the current atoms will be shifted."
                    )

parser.add_argument("--origin",
                    action="store_true",
                    help="Shift whole Box to the origin")

args = parser.parse_args()

# PROCESS FIELD ****************************************************************
finfo = dlhc.FieldFile(args.fin)
finfo.read_field()

# PROCESS CONFIG ***************************************************************
conf = dlhc.ConfigFile(finfo, args.cin)
conf.read_config()

# Shift all atoms so that center of geometry of the cell is in the origin
if args.origin:
    a = ml.vs_mult(conf.config.box.mbox[0], -0.5)
    b = ml.vs_mult(conf.config.box.mbox[1], -0.5)
    c = ml.vs_mult(conf.config.box.mbox[2], -0.5)
    conf.config.shift(a, b, c)

# Shift all atoms by the vector (defined by vectors)
if args.vectors:
    # get vectors from input
    directional_vectors = chunkIt(args.vectors, 3)
    print(directional_vectors)
    conf.config.shift(*directional_vectors)

conf.config.write_config(args.cout)
