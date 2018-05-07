#!/usr/bin/env python
from __future__ import print_function
import argparse
import matplotlib.pyplot as plt
import numpy as np

__version__ = "2017-08-07"

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="plot_energy_dat.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Plot data from paramfit output which compares " +
                "quantum, fitted and non fitted gaff parameters")

parser.add_argument("energy_dat",
                    metavar="*.dat",
                    action="store",
                    help="Paramfit output file which compares fitting results.",
                    )

args = parser.parse_args()

# read data ----------------------------------------------------------------
num = []
amber_k = []
quantum = []
initial_amber_k = []

with open(args.energy_dat) as energy_dat_in:
    energy_dat_in.next()
    column_headers = energy_dat_in.next().split()
    for line in energy_dat_in:
        line = line.split()
        num.append(int(line[0]))
        amber_k.append(float(line[1]))
        quantum.append(float(line[2]))
        initial_amber_k.append(float(line[3]))

custom_fontsize = 10
mymarkersize = 3

xyfig = plt.figure()
xyfig.canvas.set_window_title("Fit energy")
plt.title("Fit energy")
plt.xlabel("Structure", fontsize=custom_fontsize)
plt.ylabel("Energy / kcal*mol-1", fontsize=custom_fontsize)
plt.xticks(np.arange(min(num), max(num)+1, 5))
plt.plot(num, amber_k, "r-<", markersize=mymarkersize, antialiased=True, label="Fit Amber")
plt.plot(num, initial_amber_k, "b-*", markersize=mymarkersize, antialiased=True, label="Initial Amber")
plt.plot(num, quantum, "k-^", markersize=mymarkersize, antialiased=True, label="Quantum")
plt.legend(bbox_to_anchor=(0., 1., 1., .1), loc="upper right",
           borderaxespad=0., frameon=True, shadow=False, numpoints=1,
           prop={'size': 8})
plt.show()
