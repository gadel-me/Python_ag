#!/usr/bin/env python

from mpi4py import MPI
import os
import copy
import shutil as sl
import re
import argparse
import math
import numpy as np
from natsort import natsorted
#import itertools as it
import scipy.stats
from lammps import lammps, PyLammps
import Transformations as cgt
import md_elements as mde
import md_box as mdb
import ag_unify_md as agum
import ag_geometry as agm
import ag_unify_log as agul
import time
import sys

"""
Rearrange the atom order from a VESTA output file according to a given pattern.
Single atoms, which have no bonding partners are deleted during the process.
"""

s, pdb = sys.argv

lmpdat = "/home/gadelmeier/Research/FORCE_FIELDS/AMBER/Molecules/CBZ_gaff2/4.1.results/Iteration-3-2_best.lmpdat"

# read output
vesta_output = agum.Unification()
vesta_output.read_pdb(pdb)

# mark atoms which are to delete
delete_atoms = []

# find single atom fragments
for mol in vesta_output.molecules:
    if len(mol) == 1:
        idx = list(mol)
        delete_atoms.append(idx[0])

# delete all single atoms
vesta_output.delete_atoms(*delete_atoms)
vesta_output.write_xyz("tmp.xyz")
