#!/usr/bin/env python
from __future__ import print_function, division
import os
#import numpy as np
#import math
#import rmsd
#from PIL import Image

#import Transformations as cgt
#import ag_vectalg as agv
#import own_vmdfunctions as ovmd
#import ag_unify_md as agum

#import atomsel   # Replaces AtomSel and atomselection
#import axes
import color
import display
#import graphics
#import imd
#import label
#import material
import molecule
import molrep
#import mouse
#import render
#import trans
#import vmdmenu
#import Label
#import Material
#import Molecule
import VMD

"""
Load all data- and dcd-files that show the aggregation of CBZ.
"""

# change the scene
VMD.evaltcl("display shadows off")
color_display = color.get_colormap("Display")
color_display["Background"] = "white"
color.set_colormap("Display", color_display)
display.set(depthcue=0)
licor_style = "Licorice 0.1 50 50"

# get all data files
agg_maindir = "/home/gadelmeier/Unbackedup/kwz/test1/"
agg_subdirs = sorted([os.path.abspath(i) for i in os.listdir(agg_maindir) if os.path.isdir(os.path.abspath(i)) is True])
agg_cycles  = set(sorted([int(i[-1]) for i in agg_subdirs]))

print(agg_maindir, agg_subdirs, agg_cycles, os.listdir(agg_maindir))

for curcycle in agg_cycles:
    # load each stage
    for subdir in agg_subdirs:
        if "sysprep_{}".format(curcycle) in subdir:
            molecule.load("lammpsdata", subdir + "/" + "sysprep_out_{}.lmpdat".format(curcycle))
            molrep.modrep(curcycle, 0, sel="all", style=licor_style, material="AOChalky")
            break

    for subdir in agg_subdirs:
        if "quench_{}".format(curcycle) in subdir:
            molecule.read(curcycle, "dcd", subdir + "/" + "quenching_{}.dcd".format(curcycle), beg=0, end=-1, waitfor=-1)
            break

    for subdir in agg_subdirs:
        if "anneal_{}".format(curcycle) in subdir:
            molecule.read(curcycle, "dcd", subdir + "/" + "annealing_{}.dcd".format(curcycle), beg=0, end=-1, waitfor=-1)
            break

    for subdir in agg_subdirs:
        if "requench_{}".format(curcycle) in subdir:
            molecule.read(curcycle, "dcd", subdir + "/" + "requench_out_{}.dcd".format(curcycle), beg=0, end=-1, waitfor=-1)
            break
