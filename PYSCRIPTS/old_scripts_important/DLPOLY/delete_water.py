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

parser.add_argument("-fin",
                    "--dlpolyfld",
                    dest="fin",
                    metavar="xy.dlpolyfld",
                    action="store",
                    required=True,
                    help="DL_POLY-FIELD-File, belonging " +
                         "to xy.dlpolycfg/xy.dlpolyhistory",
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

#parser.add_argument("-cin2",
#                    "--config_in2",
#                    dest="cin2",
#                    metavar="xy.dlpolycfg",
#                    default=None,
#                    action="store",
#                    help="DL_POLY-CONFIG- or -REVCON-File to add cell to. " +
#                    "Cell vectors must be present in File and File must be " +
#                    "unwrapped!",
#                    )
#
#parser.add_argument("-fin2",
#                    "--dlpolyfld2",
#                    dest="fin2",
#                    metavar="xy.dlpolyfld",
#                    action="store",
#                    required=True,
#                    help="DL_POLY-FIELD-File, belonging " +
#                         "to xy.dlpolycfg/xy.dlpolyhistory",
#                    )

args = parser.parse_args()

# PROCESS FIELD ****************************************************************
finfo = dlhc.FieldFile(args.fin)
finfo.read_field()

# PROCESS CONFIG ***************************************************************
conf = dlhc.ConfigFile(finfo, args.cin)
conf.read_config()

# Delete 312 water atoms
for i in xrange(312):
    del conf.config.molecule_types[1].molecules[i]

# Get water wires
# PROCESS FIELD ****************************************************************
#finfo2 = dlhc.FieldFile(args.fin2)
#finfo2.read_field()
#
# PROCESS CONFIG ***************************************************************
#conf2 = dlhc.ConfigFile(finfo2, args.cin2)
#conf2.read_config()
#
# Extract the water wires and add them to conf
#for cur_molecule in conf2.config.molecule_types[1].molecules:
#    conf.config.molecule_types[2].molecules.append(cur_molecule)
#
conf.config.refresh_natms()
print(len(conf.config.molecule_types[0].molecules))
print(len(conf.config.molecule_types[1].molecules))
print(len(conf.config.molecule_types[2].molecules))
conf.config.write_config(config_out=args.cout)
