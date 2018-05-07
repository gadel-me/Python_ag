#!/home/angad/Programs/anaconda/bin/python -i

from __future__ import print_function, division
import os
import tqdm
#import numpy as np
import sys
sys.path.append("/home/gadelmeier/Python/myPYMODULES/DLPOPLY_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/OTHER_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/MATH_MODULES/")
import DL_HISTORY_class8_new_approaches_8_4 as dlhc
import argparse

program_version = "Unfold DL_POLY_HISTORY, build 0.95"


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

# Argument parsing stuff *******************************************************
parser = argparse.ArgumentParser(prog="0.1.Unwrap_History.py",
                                 description="Rotate and/or unwrap " +
                                 "DL_POLY-CONFIG- and/or DL_POLY-HISTORY-File")

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

parser.add_argument("-hin",
                    "--dlpolyhist",
                    dest="dlpolyhst_in",
                    metavar="xy.dlpolyhistory",
                    default=None,
                    action="store",
                    help=("DL_POLY HISTORY-File with " +
                          "box vectors for each time step"),
                    )

parser.add_argument("-hout",
                    "--dlpolyhist_out",
                    dest="dlpolyhst_out",
                    metavar="out.dlpolyhistory",
                    default=None,
                    action="store",
                    help="Filename of modified xy.dlpolyhistory",
                    )

parser.add_argument("-fld",
                    "--dlpolyfld",
                    metavar="xy.dlpolyfld",
                    action="store",
                    required=True,
                    help="DL_POLY-FIELD-File, belonging " +
                         "to xy.dlpolycfg/xy.dlpolyhistory",
                    )

#TODO STATIS-TIMESTEPS == HISTORY-TIMESTEPS
parser.add_argument("-sts",
                    "--dlstatis",
                    metavar="xy.dlstatis",
                    action="store",
                    help="DL_PLOLY-STATIS-File, belonging " +
                         "to xy.dlpolycfg/xy.dlpolyhistory." +
                         "Check if the angles are maintained after box rotation. "
                    )

parser.add_argument("--rotate",
                    action="store_true",
                    help=("Rotate the whole box so it can be understood" +
                          "by VMD and LAMMPS. " +
                          "Box-Vector a == x-axis " +
                          "Box-Vector b in xy-plane " +
                          "Box-Vector c in pos. z-direction"),
                    )

parser.add_argument("--unwrap",
                    action="store_true",
                    help=("Unwrap the whole box. " +
                          "Atoms will be outside of periodic boundaries." +
                          "***Note: It is possible that Molecules " +
                          "get shifted unreasonably. If there is no reference " +
                          "definded, first frame of xy.dlpolyhistory is used" +
                          "as reference."),
                    )

parser.add_argument("-ref",
                    "--reference_frame",
                    dest="refcon",
                    metavar="timestep|'config'",
                    default=None,
                    action="store",
                    help=("Unwrap each frame of " +
                          "xy.dlpolyhistory according to given frame. " +
                          "If argument is 'config', use xy.dlpolycfg " +
                          "as reference"),
                    )

parser.add_argument("--start",
                    dest="ts_start",
                    metavar="timestep",
                    default=None,
                    type=int,
                    action="store",
                    help="Start processing xy.dlpolyhistory with this timestep",
                    )

parser.add_argument("--stop",
                    dest="ts_stop",
                    metavar="timestep",
                    default=None,
                    type=int,
                    action="store",
                    help="Stop processing xy.dlpolyhistory beyond this timestep",
                    )

parser.add_argument("--writeframe",
                    metavar="timestep",
                    default=None,
                    action="store",
                    help="Write given timestep from xy.dlpolyhist to dlpolycfg.")

parser.add_argument("--version",
                    help="Version information",
                    action="version",
                    version=program_version
                    )

args = parser.parse_args()

# Check if args are defined in a reasonable way

# PROCESS FIELD ****************************************************************
finfo = dlhc.FieldFile(args.dlpolyfld)
finfo.read_field()

# PROCESS CONFIG/REVCON ********************************************************
# skip if not interested in CONFIG
if args.dlpolycfg_in is not None:
    # ask to delete duplicate files
    delete_file(args.dlpolycfg_out)
    conf = dlhc.ConfigFile(finfo, args.dlpolycfg_in)
    conf.read_config()

    # unwrap box
    if args.unwrap:
        # mass_center_shift => push unwrap to certain direction
        conf.config.unwrap_box(finfo,
                               mass_center_shift=[10, 10, 10])

    # rotate box
    if args.rotate:
        conf.config.rotate_box_for_vmd()

    # write new CONFIG
    if args.dlpolycfg_out is not None:
        conf.config.levcfg = 0  # write only coordinates
        conf.config.write_config(config_out=args.dlpolycfg_out,
                                 file_header=conf.config_header)

# PROCESS HISTORY **************************************************************
# Quit if not interested in HISTORY or no arguments with HISTORY
if args.dlpolyhst_in is None or (args.unwrap is False and args.rotate is False):
    if not args.dlpolycfg_in:
        print("Nothing to do.\nExiting.")
        print()
        exit(100)
    else:
        exit(0)

# ask to delete duplicate files
delete_file(args.dlpolyhst_out)

# Gather information on all timesteps within the HISTORY
hist = dlhc.HistoryFile(finfo, args.dlpolyhst_in)

# Read file until ts_stop appears (read all if ts_stop is None)
hist.read_nsteps(ts_stop=args.ts_stop, debug=False)

# Get size of one timestep
timestep_size = hist.size_of_timestep()
# Use only timesteps from given start or
# start with 1st one if no ts_start was given
process_keys = hist.nstep_keys[args.ts_start:]

# Define reference frame (only if unwrapping is checked)...
if args.unwrap:

    try:
        args.refcon = int(args.refcon)
    except (TypeError, ValueError):
        pass

    # ...by timestep
    if isinstance(args.refcon, int):
        args.refcon = int(args.refcon)
        hist.read_history(ts_start=args.refcon, ts_stop=args.refcon)
        ref_conf = hist.timesteps[args.refcon].config
    # ...by config
    elif args.refcon == "config" and args.dlpolycfg_in is not None:
        ref_conf = conf.config
    # ...by default == 1st timestep
    elif args.refcon is None:
        hist.read_history(ts_start=process_keys[0], ts_stop=process_keys[0])
        ref_conf = hist.timesteps[hist.nstep_keys[0]].config
    # ...without args.dlpolycfg_in
    elif args.dlpolycfg_in is None:
        print("***Error: Please define xy.dlpolycfg for reference to work\n")
        parser.print_help()
        exit(302)
    # ...using something undefined
    else:
        print("***Error: {}:\tInvalid argument! Choose a timestep or 'config'!\n".format(
            args.refcon)
        )
        parser.print_help()
        exit(303)

# better than re-reading file after every instance of 'HistoryFile' all over again
opened_history = open(args.dlpolyhst_in)

# stuff for making durance visible
pbar = tqdm.tqdm(process_keys)

# large history files must be read this way, else the memory-overhead is
# too massive
for istep, tqdmstep in zip(process_keys, pbar):

    # creating new instances of the class "HistoryFile" makes it possible
    # to fit large files in memory
    cur_hstep = dlhc.HistoryFile(finfo, opened_history)
    # only one frame at time
    cur_hstep.read_history(ts_start=istep, ts_stop=istep)

    # first unwrap each box according to reference...
    if args.unwrap:
        cur_hstep.timesteps[istep].config.unwrap_box(
            finfo)
        # ...then compare each frame to the unfold box (so there is no
        # unwrapping bias needed)
        # "max_discr" is the maximal allowed distance between the two
        # center of masses that are compared currently; makes the comparison run
        # faster because not every molecule has to be tested against each other
        cur_hstep.timesteps[istep].config.shift_all_mols_to_ref(
            ref_conf,
            max_discr=10)

    # rotate box so vmd can understand it
    if args.rotate:
        cur_hstep.timesteps[istep].config.rotate_box_for_vmd()

    # write (modified) frame to config
    if (args.writeframe == istep) or (args.writeframe == "all"):
        cur_config_out = "{:s}step_{:d}.dlpolyconfig".format(
            args.dlpolyhst_in,
            istep
        )
        cur_hstep.timesteps[istep].config.write_config(
            config_out=cur_config_out,
            file_header=cur_config_out,
            debug=False
        )

    # write modified history
    if args.dlpolyhst_out is not None:
        # every step is appended
        cur_hstep.write_history(history_out=args.dlpolyhst_out)

    # correct current position of file:
    # method "read_history" does not stop until having read the following
    # "TIMESTEP"-line -> file must be rewound to a point before that line

    # get current position
    cur_hist_pos = opened_history.tell()
    # cur_hist_pos must be large enough cause file-read is buffered
    cur_hist_pos -= timestep_size
    # rewind by cur_hist_pos (to ensure to get before the "TIMESTEP"-line)
    opened_history.seek(cur_hist_pos)
    # redefine the reference (to be the current config/current frame)
    ref_conf = cur_hstep.timesteps[istep].config
    # stuff for making durance visible
    pbar.set_description("Processing %s" % istep)
    pbar.update(1)

opened_history.close()
