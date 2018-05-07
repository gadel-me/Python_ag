#!/usr/bin/env python
from __future__ import print_function
import argparse
import os
import subprocess32 as sp32
import time

"""
Runs the script "gau_run.py" for each gaussian input file in the given folder.
Intended to start many gaussian jobs on a slurm cluster for parametrization/ single
point energy calculation.
"""


# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="gau_sbatch.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Script to start subscript for executing gaussian calculations.")

parser.add_argument("gau_folder",
                    help="Folder in which the gaussian input-files (*.gau) " +
                         "are stored.")

parser.add_argument("-gau_run",
                    metavar="gau_run.py",
                    help="Script to start the actual calculation")

args = parser.parse_args()

# read input-files
gau_input_files = os.listdir(args.gau_folder)
# only use input files
gau_input_files = [i for i in gau_input_files if i.endswith(".gau")]

# submit all jobs
for gau_in in gau_input_files:
    # get full path of input file
    gau_in = "{}/{}".format(args.gau_folder, gau_in)

    sp32.call(["sbatch", "--job-name", gau_in.lstrip("./"),
               args.gau_run, gau_in], bufsize=-1)

    # debugging
    #python = "/home/gadelmeier/anaconda2/bin/python"
    #sp32.call([python, args.gau_run, gau_in], bufsize=-1)
    #break
    time.sleep(0.1)
