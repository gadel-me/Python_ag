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
#import my_linalg3 as ml
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

conf2 = dlhc.ConfigFile(finfo, args.cin)
conf2.read_config()

del conf2.config.molecule_types[0].molecules
conf2.config.molecule_types[0].molecules = []

wrong_atoms = [145, 146, 757, 888,
               1250, 4128, 5489, 7630,
               7632, 9456, 9685, 9710,
               10073, 10074, 10608, 11644,
               10193, 8174, 3227, 12654,
               7693, 7695, 10367, 12714,
               12222, 11167, 4595, 7543,
               7545, 11086, 14070]
molecule_ids = []

print(len(conf.config.molecule_types[0].molecules))

for cur_mol in conf.config.molecule_types[0].molecules:
    keep_molecule = True
    for cur_atom in cur_mol.atoms:

        if cur_atom.index in wrong_atoms:
            keep_molecule = False

    if keep_molecule:
        conf2.config.molecule_types[0].molecules.append(cur_mol)

print(len(conf2.config.molecule_types[0].molecules))
conf2.config.write_config(args.cout)
