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

args = parser.parse_args()

# PROCESS FIELD ****************************************************************
finfo = dlhc.FieldFile(args.fin)
finfo.read_field()

# PROCESS CONFIG ***************************************************************
conf = dlhc.ConfigFile(finfo, args.cin)
conf.read_config()

# Half of height of H2O-Cell which was expanded by factor 1.05
h2o_1_2 = "0.000000061286      0.000000106150     34.144500500000".split()
h2o_1_2 = [float(i)*1.08 for i in h2o_1_2]
cbzii = "0.000000122572      0.000000212300     68.289001000000".split()
cbzii = [float(i) for i in cbzii]

# Everything has to be shifted by the expanded H2O cell
direction = h2o_1_2
Tm = ml.translation_matrix(direction)

for cur_moltype in conf.config.molecule_types:
    for cur_molecule in cur_moltype.molecules:
        for cur_atom in cur_molecule.atoms:
            cur_atom.xxxyyyzzz = ml.mv_mult(Tm, cur_atom.xxxyyyzzz)

# The right box height has to be defined
right_height = ml.add_vectors(h2o_1_2, ml.vs_mult(h2o_1_2, 1.07), cbzii)
conf.config.box.mbox[2] = right_height

conf.config.write_config(args.cout)
