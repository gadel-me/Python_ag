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
#import copy


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
parser = argparse.ArgumentParser(prog="cut_cell_from_cell.py",
                                 description="Cut a new cell from a given " +
                                             "cell in DLPOLY-CONFIG-Format.")

# create some mutually exclusive groups
group_a = parser.add_mutually_exclusive_group(required=True)
group_b = parser.add_mutually_exclusive_group(required=True)
group_c = parser.add_mutually_exclusive_group(required=True)

parser.add_argument("-cin",
                    "--config_in",
                    dest="cin",
                    metavar="xy.dlpolycfg",
                    default=None,
                    action="store",
                    help="DL_POLY-CONFIG- or -REVCON-File to cut cell from. " +
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

parser.add_argument("-fld",
                    "--dlpolyfld",
                    dest="fld",
                    metavar="xy.dlpolyfld",
                    action="store",
                    required=True,
                    help="DL_POLY-FIELD-File, belonging " +
                         "to xy.dlpolycfg/xy.dlpolyhistory",
                    )

group_a.add_argument("-ca",
                     "--cartesian_a",
                     dest="ca",
                     nargs=3,
                     metavar=("ax", "ay", "az"),
                     type=float,
                     default=None,
                     action="store",
                     help="Cell vector 'a' of cell which is cut from " +
                          "box in Cartesian coordinates."
                     )

group_b.add_argument("-cb",
                     "--cartesian_b",
                     dest="cb",
                     nargs=3,
                     metavar=("bx", "by", "bz"),
                     type=float,
                     default=None,
                     action="store",
                     help="Cell vector 'b' of cell which is cut from " +
                          "box in Cartesian coordinates."
                     )

group_c.add_argument("-cc",
                     "--cartesian_c",
                     dest="cc",
                     nargs=3,
                     metavar=("cx", "cy", "cz"),
                     type=float,
                     default=None,
                     action="store",
                     help="Cell vector 'c' of cell which is cut from " +
                          "box in Cartesian coordinates."
                     )

group_a.add_argument("-a",
                     "--cell_vector_a",
                     dest="a",
                     metavar=("Cell Vector a"),
                     type=float,
                     default=None,
                     action="store",
                     help="Cell vector 'a' of cell which is cut from " +
                          "box."
                     )

group_b.add_argument("-b",
                     "--cell_vector_b",
                     dest="b",
                     metavar=("Cell Vector b"),
                     type=float,
                     default=None,
                     action="store",
                     help="Cell vector 'b' of cell which is cut from " +
                          "box."
                     )

group_c.add_argument("-c",
                     "--cell_vector_c",
                     dest="c",
                     metavar=("Cell Vector c"),
                     type=float,
                     default=None,
                     action="store",
                     help="Cell vector 'c' of cell which is cut from " +
                          "box."
                     )

parser.add_argument("-alpha",
                    "--cell_angle_alpha",
                    dest="alpha",
                    metavar="Cell Angle alpha",
                    type=float,
                    default=None,
                    action="store",
                    help="Cell angle 'alpha' of cell which is cut from " +
                         "box in Cartesian coordinates in radians."
                    )

parser.add_argument("-beta",
                    "--cell_angle_beta",
                    dest="beta",
                    metavar="Cell Angle beta",
                    type=float,
                    default=None,
                    action="store",
                    help="Cell angle 'beta' of cell which is cut from " +
                         "box in Cartesian coordinates in radians."
                    )

parser.add_argument("-gamma",
                    "--cell_angle_gamma",
                    dest="gamma",
                    metavar="Cell Angle gamma",
                    type=float,
                    default=None,
                    action="store",
                    help="Cell angle 'gamma' of cell which is cut from " +
                         "box in Cartesian coordinates in radians."
                    )

args = parser.parse_args()

# ask to delete duplicate files
if args.cout is not None:
    delete_file(args.cout)

# PROCESS FIELD ****************************************************************
finfo = dlhc.FieldFile(args.fld)
finfo.read_field()

# PROCESS CONFIG ***************************************************************
conf = dlhc.ConfigFile(finfo, args.cin)
conf.read_config()

# CUT CONFIG *******************************************************************
# Create a new config-instance that stores only those molecules that are cut
# from the main config-instance

# Create new box vectors
if args.a or args.b or args.c:
    (a_cartesian,
     b_cartesian,
     c_cartesian,
     M_fract2cart) = ml.vt_cell2box(args.a, args.b, args.c,
                                    args.alpha, args.beta, args.gamma)
else:
    a_cartesian = args.ca
    b_cartesian = args.cb
    c_cartesian = args.cc

cut_conf = dlhc.ConfigFile(finfo, args.cin)
cut_conf.config_header = conf.config_header + " cut"
cut_conf.field_info = conf.field_info
cut_conf.dl_config = conf.dl_config
cut_conf.config.molecule_types = finfo.molecule_types[:]

# Box stuff
cut_conf.config.imcon = conf.config.imcon
cut_conf.config.levcfg = conf.config.levcfg
cut_conf.config.box = dlhc.Box()
cut_conf.config.box.mbox[0] = a_cartesian
cut_conf.config.box.mbox[1] = b_cartesian
cut_conf.config.box.mbox[2] = c_cartesian
print(a_cartesian, b_cartesian, c_cartesian)
# STUFF FOR CARTESIAN TO FRACTIONAL TRANSFORMATION *****************************

# Generate cell angles if no angles were given
if args.alpha is None and not(args.a or args.b or args.c):
    args.alpha = ml.get_angle_between(args.cb, args.cc)
if args.beta is None and not(args.a or args.b or args.c):
    args.beta = ml.get_angle_between(args.ca, args.cc)
if args.gamma is None and not(args.a or args.b or args.c):
    args.gamma = ml.get_angle_between(args.ca, args.cb)

# Convert Cartesian box vectors to cell vectors
if args.ca:
    args.a = ml.get_magnitude(args.ca)
if args.cb:
    args.b = ml.get_magnitude(args.cb)
if args.cc:
    args.c = ml.get_magnitude(args.cc)

# Get fractional matrix
M_fract = ml.M_cart2fract(args.a, args.b, args.c,
                          args.alpha, args.beta, args.gamma)

# Start re-indexing the atom numbers
new_atom_index = 0

for conf_moltype, cut_moltype in zip(conf.config.molecule_types,
                                     cut_conf.config.molecule_types):
    for cur_mol in conf_moltype.molecules:
        delete_molecule = False

        for cur_atom in cur_mol.atoms:
            # Transform the Cartesian coordinates to fractional coordinates
            # of the new cell
            cur_atom_fract_coords = ml.coords_cart2fractional(M_fract, cur_atom.xxxyyyzzz)

            # Fract. coords. within the cell are between 0 and 1 for
            # each cell vector; e.g. 0.5*a, 0.3*b, 0.02*c
            if not (cur_atom_fract_coords[0] >= 0 and
                    cur_atom_fract_coords[0] <= 1) or\
               not (cur_atom_fract_coords[1] >= 0 and
                    cur_atom_fract_coords[1] <= 1) or\
               not (cur_atom_fract_coords[2] >= 0 and
                    cur_atom_fract_coords[2] <= 1):
                delete_molecule = True
                break

        if not delete_molecule:
            #print("Keeping molecule id: ", cur_mol.mol_id)
            cut_moltype.molecules.append(cur_mol)

#TODO Write corresponding field file (just change number of atoms in finfo)
# Refresh information
cut_conf.config.refresh_natms()
cut_conf.config.levcfg = 0
cut_conf.config.write_config(config_out=args.cout)
