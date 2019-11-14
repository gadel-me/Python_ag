#!/usr/bin/env python

import os
import numpy as np
import math
import rmsd
from PIL import Image

import Transformations as cgt
import ag_vectalg as agv
import own_vmdfunctions as ovmd
import ag_unify_md as agum

import atomsel   # Replaces AtomSel and atomselection
import axes
import color
import display
import graphics
import imd
import label
import material
import molecule
import molrep
import mouse
import render
import trans
import vmdmenu
import Label
import Material
import Molecule
import VMD

"""
Load cbz-ab-initio, cbz-gaff and cbz-gaff-opt average-structure-files,
change their representation and take a snapshot of the scene.
"""

home = os.environ["HOME"]

# structure files
cbz_gOpt         = "{}/Research/FORCE_FIELDS/AMBER/Molecules/CBZ_gaff2/1.resp_charges/2-geomOpt/CBZ_gau_Opt_log.xyz".format(home)
#cbz_gaff2_lmpdat = "{}/Research/FORCE_FIELDS/AMBER/Molecules/CBZ_gaff2/4.2.test_MDs/Iteration_3-2/Iteration-3-2_best.lmpdat".format(home)
#cbz_gaff2_dcd    = "{}/Research/FORCE_FIELDS/AMBER/Molecules/CBZ_gaff2/4.2.test_MDs/Iteration_3-2/Iteration-3-2_best_in_vacuo.dcd".format(home)
#cbz_gaff_lmpdat  = "{}/Research/FORCE_FIELDS/AMBER/Molecules/CBZ_gaff2/4.2.test_MDs/Standard_gaff/CBZ.lmpdat".format(home)
#cbz_gaff_dcd     = "{}/Research/FORCE_FIELDS/AMBER/Molecules/CBZ_gaff2/4.2.test_MDs/Standard_gaff/Standard_gaff_in_vacuo.dcd".format(home)
cbz_gaff2_lmpdat = "{}/Research/PAPERS/2017_12/CBZ_FF-evaluation/results/average_structures/cbz_gaff2_average.lmpdat".format(home)
cbz_gaff_lmpdat  = "{}/Research/PAPERS/2017_12/CBZ_FF-evaluation/results/average_structures/cbz_gaff_average.lmpdat".format(home)

# change the scene
VMD.evaltcl("display shadows off")
color_display = color.get_colormap("Display")
color_display["Background"] = "white"
color.set_colormap("Display", color_display)
display.set(depthcue=0)

# load mol and change its representation
# geometry optimized cbz
licor_style = "Licorice 0.1 50 50"
sel_no_h    = "not element H"

if molecule.exists(0) != 1:
    molecule.load("xyz", cbz_gOpt)
    molrep.modrep(0, 0, sel=sel_no_h, style=licor_style, material="AOChalky", color="ColorID 7")

if molecule.exists(1) != 1:
    #molecule.load("lammpsdata", cbz_gaff2_lmpdat, "dcd", cbz_gaff2_dcd)
    molecule.load("lammpsdata", cbz_gaff2_lmpdat)
    molrep.modrep(1, 0, sel=sel_no_h, style=licor_style, material="AOChalky", color="ColorID 0")

if molecule.exists(2) != 1:
    #molecule.load("lammpsdata", cbz_gaff_lmpdat, "dcd", cbz_gaff_dcd)
    molecule.load("lammpsdata", cbz_gaff_lmpdat)
    molrep.modrep(2, 0, sel=sel_no_h, style=licor_style, material="AOChalky", color="ColorID 1")

# calculate the average structure of each molecule


# draw half circle for angle representation
ovmd.draw_angle(0, 0, (15, 24, 7), canvas=True, resolution=50, color="red")
ovmd.draw_angle(0, 0, (15, 24, 7), canvas=False, resolution=50, color="blue")
ovmd.draw_angle(1, 0, (15, 24, 7), canvas=True, resolution=50, color="red")
ovmd.draw_angle(1, 0, (15, 24, 7), canvas=False, resolution=50, color="blue")
VMD.evaltcl("display resetview")

# rotate camera
trans.rotate("z", 180)
trans.rotate("x", -35)

# zoom out
for _ in range(8):
    trans.scale(0.833)

# rendering
render_scene = True  # just for testing stuff
if render_scene is True:
    disp_default_size = display.get("size")
    display.set(size=[4000, 4000])  # higher resolution
    image_out = "{}/Research/PAPERS/2017_12/CBZ_FF-evaluation/vmd_scenes/cbz_comparison".format(home)
    render.render("TachyonLOptiXInternal", "{}.ppm".format(image_out))
    im = Image.open("{}.ppm".format(image_out))
    im.save("{}.png".format(image_out), format="PNG")
    display.set(size=disp_default_size)
    os.remove("{}.ppm".format(image_out))
