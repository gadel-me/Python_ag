#!/usr/bin/env python
from __future__ import print_function, division
import math
import numpy as np
import os
import shutil as sl
import argparse
from mpi4py import MPI
from lammps import lammps, PyLammps
import ag_unify_md as agum
import ag_geometry as agm
import Transformations as cgt

# MPI initialization -----------------------------------------------------------
comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="kawska_zahn.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Kawska-Zahn-Approach for accelerated crystallization simulations.")

# mutually exclusive groups
starters = parser.add_mutually_exclusive_group(required=True)

starters.add_argument("-lmpdat",
                      metavar="*.lmpdat",
                      default=None,
                      help="lammps' data-file"
                      )

parser.add_argument("-mols2dock",
                    nargs="*",
                    metavar="*.lmpdat",
                    help="lammps' data-file with atom-types and coordinates \
                         for one single molecule to add. Atom types have to be \
                         defined already by the first data/restart file loaded!")

starters.add_argument("-lmprst",
                      metavar="*.lmprst",
                      default=None,
                      help="lammps' restart file",
                      )

parser.add_argument("-solventmols",
                    nargs="*",
                    metavar="*.lmpdat",
                    default=None,
                    help="Create a solvent box for MD-Simulation. \
                          Data-file with one single solvent molecule to add."
                    )

parser.add_argument("-set",
                    metavar="*.lmpcfg",
                    required=True,
                    help="lammps' input-file/-script with basic simulation " +
                         "settings"
                    )

parser.add_argument("-steps",
                    type=int,
                    default=5000,
                    help="number of timesteps for each run"
                    )

parser.add_argument("-logsteps",
                    type=int,
                    default=1000,
                    help="log thermodynamic-, steps- and restart-files every" +
                         "'logsteps' steps"
                    )

parser.add_argument("-log",
                    default=None,
                    help="file that logs the current state of the \
                          kawska-zahn approach (needed for restarting stuff)")

parser.add_argument("-gpu",
                    default=False,
                    action="store_true",
                    help="utilize lammps' GPU package.",
                    )

args = parser.parse_args()

# LAMMPS BASIC -----------------------------------------------------------------
lmp = lammps()  # (re-)create lammps-instance
LMP = PyLammps(ptr=lmp)

# load gpu package
if args.gpu:
    # neighbor list building on gpu not available for triclinic boxes (neigh no)
    lmp.command("package gpu 1 neigh no")
    lmp.command("suffix gpu")


# load general settings (units, boundary, dimension, atom_style, etc. pp.)
lmp.file(args.set)
lmp.command("box tilt large")  # do not care about too tilted boxes

# load the restart or data file (prefer restart-file over data-file)
if args.lmprst:
    lmp.command("read_restart {}".format(args.lmprst))
elif args.lmpdat:
    lmp.command("read_data {}".format(args.lmpdat))
else:
    raise RuntimeError("Neither restart nor data file defined!")

# define thermo output
Thermargs = ["step", "temp", "press", "vol", "density",
             "cella", "cellb", "cellc", "cellalpha", "cellbeta", "cellgamma",
             "atoms", "bonds", "angles",
             "etotal", "pe", "evdwl", "ecoul", "ebond", "eangle", "edihed", "eimp",
             "enthalpy"]
lmp.command("thermo_style custom " + " ".join(Thermargs))
lmp.command("thermo_modify lost warn flush yes")
lmp.command("thermo {}".format(int(math.ceil(args.logsteps/2))))

# docking ----------------------------------------------------------------------

# new coordinates for molecule from data file to be placed around given atoms
if rank == 0:
    mols2add = []
    for cmol2dock in args.mols2dock:
        cur_dockmol = agum.Unification()
        cur_dockmol.read_data(cmol2dock)

        # shift molecue to origin
        cur_dockmol_idxs = xrange(len(cur_dockmol.atoms))
        com_dockmol = cur_dockmol.get_com(0, cur_dockmol_idxs)

        if com_dockmol != [0, 0, 0]:
            Tcom = cgt.translation_matrix(com_dockmol)
            cur_dockmol.mm_atm_coords(0, Tcom, False, cur_dockmol_idxs)

        # molecule rotation (define internal molecule coord system before)
        # internal coord sys
        com_bzrght = cur_dockmol.get_com(0, 2, 4, 6, 7, 8)
        com_bzlft  = cur_dockmol.get_com(14, 15, 16, 18, 20, 22)
        vt_bzs = np.array(com_bzlft) - np.array(com_bzrght)
        vt_nc  = cur_dockmol.ts_coords[0][26] - cur_dockmol.ts_coords[0][25]
        bzx, bzy, bzz = agm.get_coord_sys(vt_bzs, vt_nc)
        # rotation about arbitrary vector
        Marb_rotx = agm.arb_rot_matrix(bzx)
        Marb_roty = agm.arb_rot_matrix(bzy)
        Marb_rotz = agm.arb_rot_matrix(bzz)
        Mrot_xy = np.matmul(Marb_rotx, Marb_roty)
        Mrot = np.matmul(Marb_rotz, Mrot_xy)

        # shift molecule to random point on sphere
        rn_pos = agum.points_on_sphere(npoints=1, radius=10)
        Tm = cgt.translation_matrix(rn_pos)
        cur_dockmol.mm_atm_coords(0, Tm, False, cur_dockmol_idxs)
        # 'save' molecule
        mols2add.append(cur_dockmol)
else:
    mols2add = None

comm.bcast(mols2add, root=0)
