#!/usr/bin/env python
#SBATCH --exclusive

from __future__ import print_function
import argparse
import os
import re
import subprocess32 as sp32
import shutil as sl
import psutil
import math
import ag_gaussian as aggau


# Helper functions for cpu and memory capacities on current node ---------------
def get_real_cores():
    ht_info = os.popen('lscpu').readlines()[6]
    socket_info = os.popen('lscpu').readlines()[7]
    cores = int(re.search(r'\d+', ht_info).group(0))
    sockets = int(re.search(r'\d+', socket_info).group(0))
    return (cores*sockets)


def get_real_memory():
    mem = psutil.virtual_memory()
    mem_in_gb = int(math.floor(mem.available*9.31323e-10))
    return mem_in_gb


# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="gau_run.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Submit a gaussian calculation.")

parser.add_argument("gau_in",
                    help="Gaussian input file (.gau).")

args = parser.parse_args()

# change cpu and memory entries to what we have
gau = aggau.Gaussian(nproc=get_real_cores(), mem=get_real_memory())

# read input file
gau.read_gau(args.gau_in)

# write new output file
gau.write_gau(args.gau_in+".tmp", 0)
sl.move(args.gau_in+".tmp", args.gau_in)

# execute gaussian
g09 = "/usr/bin/g09"
sp32.call([g09, args.gau_in], bufsize=-1)

# rename gaussian-out (.gau.out -> .log)
g_out = args.gau_in + ".out"
g_log = g_out.replace(".gau.out", ".log")
sl.move(g_out, g_log)
