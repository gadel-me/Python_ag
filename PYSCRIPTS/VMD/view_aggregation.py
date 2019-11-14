#!/usr/bin/env python

import pdb
import os
import re
import pathlib

import ag_vmd

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


def get_finished_cycles(maindir):
    """
    """
    fc = []
    folders = ["{}/{}".format(maindir, i) for i in os.listdir(maindir) if os.path.isdir(i)]
    # get last cycle from directory
    for folder in folders:
        # skip all folders where anything went wrong
        if "fail" in folder:
            continue
        #print(folder)
        cycle = re.match(r'.*?([0-9]+)$', folder).group(1)
        cycle = int(cycle)
        # avoid duplicates
        if cycle not in fc:
            fc.append(cycle)

    return fc


def get_files(maindir, filepattern):
    files = []
    #pdb.set_trace()

    for cfile in pathlib.Path(maindir).glob("**/{}".format(filepattern)):
        files.append(cfile)

    str_filenames = []

    for filename in files:
        full_path = filename.resolve()
        str_filename = full_path.as_posix()
        str_filenames.append(str_filename)

    return str_filenames


if __name__ == "__main__":
    MAINDIR = "."
    ITERATIONS = None

    SELECTIONS = [
        {"style": "Lines 1.0", "color": "Name", "selection": "type 'c'", "material": "Basic1Pantone"},
        {"style": "Lines 2.0", "color": "Name", "selection": "name C15 or type hn or type o or name N2", "material": "Basic1Pantone"}
    ]

    SYSPREP_LMPDAT_RAW = "{0}/sysprep_{1}/sysprep_out_{1}.lmpdat"
    QUENCH_DCD_RAW = "{0}/quench_{1}/quench_{1}.dcd"
    EQUIL_ANNEAL_DCD_RAW = "{0}/anneal_{1}/equil_anneal_{1}.dcd"
    PROD_ANNEAL_DCD_RAW = r"[0-9]+_anneal_{0}.dcd"
    #REQUENCH_LMPDAT = "{0}/requench_{1}/requench_out_{1}.lmpdat"
    REQUENCH_DCD_RAW = "{0}/requench_{1}/requench_{1}.dcd"

    # find all iterations
    if ITERATIONS is None:
        FINISHED_CYCLES = get_finished_cycles(MAINDIR)
    else:
        FINISHED_CYCLES = ITERATIONS

    for CURCYCLE in FINISHED_CYCLES:
        # sysprep
        SYSPREP_LMPDAT = SYSPREP_LMPDAT_RAW.format(MAINDIR, CURCYCLE)
        ag_vmd.vmd_load_molecule(SYSPREP_LMPDAT, style="Lines 1.000", molid=CURCYCLE)
        molrep.set_visible(CURCYCLE, 0, False)

        # add representations
        for CSEL in SELECTIONS:
            molrep.addrep(
                CURCYCLE,
                style=CSEL["style"],
                color=CSEL["color"],
                selection=CSEL["selection"],
                material=CSEL["material"])

        # quenching
        QUENCH_DCD = QUENCH_DCD_RAW.format(MAINDIR, CURCYCLE)
        molecule.read(CURCYCLE, "dcd", QUENCH_DCD, beg=0, end=-1, waitfor=-1)

        # annealing - heating up
        EQUIL_ANNEAL_DCD = EQUIL_ANNEAL_DCD_RAW.format(MAINDIR, CURCYCLE)
        molecule.read(CURCYCLE, "dcd", EQUIL_ANNEAL_DCD, beg=0, end=-1, waitfor=-1)

        # annealing - productive
        ANNEAL_DIR = "{}/anneal_{}".format(MAINDIR, CURCYCLE)
        PROD_ANNEAL_DCDS = get_files(
            "{0}/anneal_{1}/".format(ANNEAL_DIR, CURCYCLE),
            PROD_ANNEAL_DCD_RAW.format(CURCYCLE))

        for CUR_DCD in PROD_ANNEAL_DCDS:
            print(CUR_DCD)
            molecule.read(CURCYCLE, "dcd", CUR_DCD, beg=0, end=-1, waitfor=-1)

        pdb.set_trace()

        # requenching
        REQUENCH_DCD = REQUENCH_DCD_RAW.format(MAINDIR, CURCYCLE)
        molecule.read(CURCYCLE, "dcd", REQUENCH_DCD, beg=0, end=-1, waitfor=-1)
        VMD.evaltcl("mol off {}".format(CURCYCLE))
