#!/bin/env python

from __future__ import print_function, division
#import numpy as np
import sys
import argparse
import os

# user-modules******************************************************************
home = os.environ["HOME"]
sys.path.append(home + "/Python/myPYMODULES" + "/DLPOPLY_MODULES")
sys.path.append(home + "/Python/myPYMODULES" + "/OTHER_MODULES")
sys.path.append(home + "/Python/myPYMODULES" + "/MATH_MODULES")
import DL_HISTORY_class8_new_approaches_8_4 as dlhc

# Stuff one may lookup for further reading
# remove items during loop from list
# From: http://stackoverflow.com/questions/6022764/\
# python-removing-list-element-while-iterating-over-list


# Argument parsing stuff *******************************************************
parser = argparse.ArgumentParser(prog="remove_half.py",
                                 description="Cut part of the box")

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

parser.add_argument("-z_min",
                    action="store",
                    default=None,
                    type=float,
                    help="minimal z-position to cut")

parser.add_argument("-z_max",
                    action="store",
                    default=None,
                    type=float,
                    help="maximal z-position to cut")


args = parser.parse_args()

# PROCESS FIELD ****************************************************************
finfo = dlhc.FieldFile(args.dlpolyfld)
finfo.read_field()


# PROCESS CONFIG/REVCON ********************************************************
conf = dlhc.ConfigFile(finfo, args.dlpolycfg_in)
conf.read_config()

mol_to_del = []

for cur_moltype in conf.config.molecule_types:
    if cur_moltype.molname == "TIP3P water":

        for cur_molecule in list(cur_moltype.molecules):
            delete_molecule = False
            for cur_atom in cur_molecule.atoms:

                if (cur_atom.xxxyyyzzz[2] > args.z_min) and (cur_atom.xxxyyyzzz[2] < args.z_max):
                    delete_molecule = True

            if delete_molecule:
                cur_moltype.molecules.remove(cur_molecule)

conf.config.refresh_natms()
#conf.config.rewrite_field(conf.field_info)
conf.config.write_config(config_out=args.dlpolycfg_out,
                         file_header=conf.config_header)
