#!/usr/bin/env python
#SBATCH --job-name=
#SBATCH --exclusive
#SBATCH --output=slurm.%N.%J.out
#SBATCH --error=slurm.%N.%J.err
#SBATCH --partition=24h
#SBATCH --exclude=znode[401,402,403]

from __future__ import print_function
import sys
import os
import re
import subprocess32 as sp32
import shutil as sl
import psutil
import math


def get_real_cores():
    ht_info = os.popen('lscpu').readlines()[6]
    socket_info = os.popen('lscpu').readlines()[7]
    cores = int(re.search(r'\d+', ht_info).group(0))
    sockets = int(re.search(r'\d+', socket_info).group(0))
    return str(cores*sockets)


def get_real_memory():
    mem = psutil.virtual_memory()
    mem_in_gb = int(math.floor(mem.total*9.31323e-10))
    return mem_in_gb

script, gau_in = sys.argv
g09 = "/usr/bin/g09"
gau_in_temp = gau_in + ".temp"
cores_avail = get_real_cores()
mem_avail = get_real_memory()

# Rewrite gaussian input (adjust number of cores and memory)
with open(gau_in, "r") as g_in, open(gau_in_temp, "w") as g_temp:
    for line in g_in:

        # change number of utilized cores
        if "%nproc" in line:
            line = "%nproc={}\n".format(cores_avail)
        # change number of available memory
        if "%mem" in line:
            line = "%mem={}GB\n".format(mem_avail)

        g_temp.write(line)

sl.move(gau_in_temp, gau_in)
sp32.call([g09, gau_in], bufsize=-1)

# rename gaussian-out (xyz.gau.out) to xyz.log
g_out = gau_in + ".out"
g_log = g_out.replace(".gau.out", ".log")
sl.move(g_out, g_log)
