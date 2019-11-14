#!/usr/bin/env python

import os
import argparse
import subprocess32 as sp32
import time

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="jaguar_batch_processing.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Run all jaguar input scripts."
)

parser.add_argument("-jag_folders",
                    nargs="*",
                    action="store",
                    required=True
                    )

parser.add_argument("-estimated_duration",
                    type=int,
                    default=3600,
                    help="Estimated time one run needs (in seconds)")

args = parser.parse_args()

# fetch some dirs --------------------------------------------------------------
workdir = os.getcwd()
schrodinger = os.getenv("SCHRODINGER")

for jag_folder in args.jag_folders:
    os.chdir(jag_folder)
    print("Entering directory: ", jag_folder)
    jaguar_in = "{0}.in".format(jag_folder)
    jaguar_run = sp32.call(["{}/jaguar".format(schrodinger), "run", "-TPP", "4", jaguar_in])
    os.chdir(workdir)
    print("Returning to ", workdir)
    print("Waiting for {} h (hopefully job finishes earlier)!".format(args.estimated_duration/3600))
    time.sleep(args.estimated_duration)
