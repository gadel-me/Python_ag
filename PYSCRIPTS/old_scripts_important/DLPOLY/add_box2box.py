#!/home/angad/Programs/anaconda/bin/python -i

from __future__ import print_function, division
import os
#import tqdm
#import numpy as np
import sys
sys.path.append("/home/gadelmeier/Python/myPYMODULES/DLPOPLY_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/OTHER_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/MATH_MODULES/")
import DL_HISTORY_class8_new_approaches_8_4 as dlhc
import my_linalg3 as ml
import argparse
import copy


def delete_file(filename):
    """
    Look up, if a file is already existent and ask to delete if it is.
    """
    # end function if no filename is given
    if filename is None:
        return None
    asking = True
    if os.path.isfile(filename):
        while asking:
            deletion = raw_input("***Warning: '{:s}' is already existing! ".format(filename) +
                                 "Delete file <Y/n>?\n> "
                                 )

            if deletion == "Y":
                os.remove(filename)
                asking = False
            elif deletion == "n":
                print("***Note: Rename or remove file or change name '{:s}' for this program to work.".format(
                    filename)
                )
                print("Exiting.")
                exit(001)
            else:
                print("***{:s} not a valid answer! Type 'Y' for (Y)es or 'n' for (n)o!".format(
                    deletion)
                )

program_version = "Cut CONFIG-cell, build 0.01"

# ARGUMENT PARSING STUFF *******************************************************
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

parser.add_argument("-cadd",
                    "--add_config_cell",
                    dest="cadd",
                    metavar="xy.dlpolycfg",
                    default=None,
                    action="store",
                    help="DL_POLY-CONFIG- or -REVCON-File which is added. " +
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

parser.add_argument("-fadd",
                    "--fadd",
                    dest="fld_add",
                    metavar="xy.dlpolyfld",
                    action="store",
                    required=True,
                    help="DL_POLY-FIELD-File, belonging " +
                         "to xy.dlpolycfg/xy.dlpolyhistory",
                    )

parser.add_argument("-vt",
                    dest="vt",
                    default=None,
                    action="store",
                    help="Cell vector in which direction the box will be" +
                         "added (in Cartesian coordinates)"
                    )

parser.add_argument("-buffer",
                    type=float,
                    default=None,
                    action="store",
                    help=("Multiply vt by the buffer to prevent clashes.")
                    )

args = parser.parse_args()

# ask to delete duplicate files
if args.cout is not None:
    delete_file(args.cout)

# PROCESS FIELD ****************************************************************
finfo = dlhc.FieldFile(args.fin)
finfo.read_field()

# PROCESS CONFIG ***************************************************************
conf = dlhc.ConfigFile(finfo, args.cin)
conf.read_config()

# PROCESS CELL TO ADD
finfo_add = dlhc.FieldFile(args.fld_add)
finfo_add.read_field()

cadd_conf = dlhc.ConfigFile(finfo_add, args.cadd)
cadd_conf.read_config()

# pos cell direction
if args.vt == "c":
    # Get enlarged cell vector c
    cell_c_add = ml.vs_mult(cadd_conf.config.box.mbox[2], args.buffer)
    # Get the difference in length between c(original) and c(enlarged)
    cell_c_space = ml.add_vectors(cell_c_add, [-1*i for i in cadd_conf.config.box.mbox[2]])
    buffered_cdir = ml.add_vectors(conf.config.box.mbox[2], cell_c_space)
    Tm1 = ml.translation_matrix(buffered_cdir)

cadd_1_conf = copy.deepcopy(cadd_conf)
for cur_moltype in cadd_1_conf.config.molecule_types:
    for cur_mol in cur_moltype.molecules:
        for cur_atom in cur_mol.atoms:
            cur_atom.xxxyyyzzz = ml.mv_mult(Tm1, cur_atom.xxxyyyzzz)
        conf.config.molecule_types[1].molecules.append(cur_mol)

# neg cell direction
buffered_cdir2 = [-1*i for i in ml.vs_mult(cadd_conf.config.box.mbox[2], args.buffer)]
Tm2 = ml.translation_matrix(buffered_cdir2)

cadd_2_conf = copy.deepcopy(cadd_conf)
for cur_moltype in cadd_2_conf.config.molecule_types:
    for cur_mol in cur_moltype.molecules:
        for cur_atom in cur_mol.atoms:
            cur_atom.xxxyyyzzz = ml.mv_mult(Tm2, cur_atom.xxxyyyzzz)
        conf.config.molecule_types[1].molecules.append(cur_mol)

conf.config.refresh_natms()

if args.vt == "a":
    total_length = ml.add_vectors(conf.config.box.mbox[0], cell_c_add, cell_c_add)
    conf.config.box.mbox[0] = total_length
elif args.vt == "b":
    total_length = ml.add_vectors(conf.config.box.mbox[1], cell_c_add, cell_c_add)
    conf.config.box.mbox[1] = total_length
elif args.vt == "c":
    total_length = ml.add_vectors(conf.config.box.mbox[2], cell_c_add, cell_c_add, [0.5*i for i in cell_c_space])
    conf.config.box.mbox[2] = total_length
else:
    print("Not a valid cell vector")

conf.config.levcfg = 0

# shift to origin
a_2, b_2, c_2 = [((i*-0.5) for i in j) for j in conf.config.box.mbox]
conf.config.shift(a_2, b_2, c_2, cell_c_add)
conf.config.write_config(config_out=args.cout)
