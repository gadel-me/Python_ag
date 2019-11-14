#!/usr/bin/env python

import os
import argparse
import sys
import subprocess32 as sp32
import shutil as sl

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="create_inputs.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Converts maestro output-files to amber coordinates and \
                gaussian-input-files."
)

parser.add_argument("jag_folders",
                    nargs="*",
                    help="Jaguar input folder."
                    )

parser.add_argument("-mdcrd_dir",
                    default=None,
                    help="Optional directory for amber coordinate file."
                    )

args = parser.parse_args()

# subscripts -------------------------------------------------------------------
python = sys.executable
mae2mdcrd_gau = "/home/gadelmeier/Python/PYSCRIPTS/MAESTRO/mae2mdcrd_gau.py"
mae_out = None

if args.mdcrd_dir is not None:
    args.mdcrd_dir = os.path.abspath(args.mdcrd_dir)

# convert files to gaussian input and amber coordinates
for jagfolder in args.jag_folders:
    # get absolute path of current folder
    jagfolder = os.path.abspath(jagfolder)

    # get files from folder
    maestro_files = os.listdir(jagfolder)

    # find maestro output file
    for cf in maestro_files:
        if cf.endswith(".01.mae"):
            name_maestro_out = cf
            break

    # absolute path to maestro file
    maestro_01_mae = "{}/{}".format(jagfolder, name_maestro_out)
    projectname = name_maestro_out.rstrip(".01.mae")

    # create new directory or skip if it already exists
    try:
        os.mkdir(projectname)
    except OSError:
        # skip if files are in folder
        if os.listdir(projectname) != []:
            print("Folder already exists. Skipping")
            continue
        else:
            pass

    os.chdir(projectname)
    sp32.call([python, mae2mdcrd_gau, "-mae", maestro_01_mae, "-out", projectname])

    if os.path.isdir(args.mdcrd_dir) is True:
        # move amber coordinate file
        mdcrd_src = "{}.mdcrd".format(projectname)
        mdcrd_dst = "{}/{}".format(args.mdcrd_dir, mdcrd_src)
        sl.move(mdcrd_src, mdcrd_dst)

    os.chdir(os.pardir)
