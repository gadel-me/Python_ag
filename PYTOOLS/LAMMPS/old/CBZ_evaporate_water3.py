#!/usr/bin/env python

from mpi4py import MPI
from lammps import lammps
import os
import sys
import argparse

home = os.getenv("HOME")
sys.path.append(home + "/Python/myPYMODULES/DLPOPLY_MODULES")
sys.path.append(home + "/Python/myPYMODULES/OTHER_MODULES")
sys.path.append(home + "/Python/myPYMODULES/MATH_MODULES")
sys.path.append(home + "/Python/myPYMODULES/LAMMPS")
import lmp_module_ga as lmga
import shutil as sl

# import collections
# import numpy as np
# import math
# import warnings
# import re
# import time
# import ctypes

# ==============================================================================#
# ARGUMENT PARSING=============================================================#
# ==============================================================================#
parser = argparse.ArgumentParser(
    prog="vaporize.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Delete molecules above " + "certain border from the simulation.",
)

parser.add_argument(
    "-fset",
    metavar="*.lmpcfg",
    action="store",
    required=True,
    help="lammps' input-file/-script with basic simulation\n"
    + "settings\n"
    + "http://lammps.sandia.gov/doc/Section_commands.html",
)

parser.add_argument(
    "-fdat",
    metavar="*.lmpdat",
    nargs="*",
    action="store",
    required=True,
    help="lammps' data-file(s);\n"
    + "box-, force field-, "
    + "topology-parameters must be included!\n"
    + "http://lammps.sandia.gov/doc/2001/data_format.html",
)

parser.add_argument(
    "-fpc",
    metavar="*.lmpcfg",
    action="store",
    required=True,
    help="lammps' input-file/-script with pair coefficients\n"
    + "http://lammps.sandia.gov/doc/pair_modify.html",
)

parser.add_argument(
    "-rst",
    "--restart_file",
    dest="rst",
    metavar="*.lmprst",
    default=None,
    action="store",
    help="lammps' restart file (load instead of data-file)\n"
    + "http://lammps.sandia.gov/doc/read_restart.html",
)

parser.add_argument(
    "-o",
    "--outputname",
    dest="o",
    default="DEFAULTNAME",
    action="store",
    help="output file names",
)

parser.add_argument(
    "-iter",
    "--iterations",
    dest="iter",
    type=int,
    default=1,
    metavar=1,
    help="number of simulation iterations",
)

parser.add_argument(
    "-steps",
    type=int,
    metavar=1000,
    default=5000,
    action="store",
    help="number of timesteps for each run",
)

parser.add_argument(
    "-logsteps",
    type=int,
    default=1000,
    help="log thermodynamic-, steps- and restart-files every" + "logsteps steps",
)

parser.add_argument(
    "-zmax",
    type=float,
    metavar=55.0,
    default=55.0,
    help="***WARING: UNDER CONSTRUCTION, MAX Z-COORD FOR ATOM.",
)

parser.add_argument(
    "-gpu", default=False, action="store_true", help="utilize lammps' GPU package.",
)

args = parser.parse_args()

# ==============================================================================#
# MPI-INITIALIZATION ==========================================================#
# ==============================================================================#
comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id within a communicator

# ==============================================================================#
# GENERAL SETTINGS ============================================================#
# ==============================================================================#
# * Trajectory-/restart-files
Dsteps = args.logsteps  # write frame every Dsteps
Rsteps = args.logsteps  # write restart every Rsteps
# * Temperature
Tstart = 300.0
Tstop = Tstart
thermsteps = args.logsteps  # write thermodynamic info every Thermsteps
# boundaries (for deletion)
zmax1 = args.zmax  # upper boundary
zmax2 = -zmax1  # lower boundary
first_run = True  # set to True if script was initialized the first time
new_run = False  # set to True if a run is restarted (atoms got deleted)

# ==============================================================================#
# ITERATION-SETTINGS ==========================================================#
# ==============================================================================#

for cur_iteration in range(args.iter):
    # (re-)create lammps-instance
    lmp = lammps()

    # ==========================================================================#
    # PREREQUESITES ===========================================================#
    # ==========================================================================#
    # load gpu package
    if args.gpu:
        # neighbor list building on gpu not available for triclinic boxes (neigh no)
        lmp.command("package gpu 1 neigh no")
        lmp.command("suffix gpu")

    # load settings
    lmp.file(args.fset)

    # ==========================================================================#
    # (RE-)LOAD DATA/RESTART-FILE =============================================#
    # ==========================================================================#
    data_index = []  # placeholder
    data_shift = ""  # placeholder

    # use files from command line if first run
    if first_run:

        if args.rst:  # prefer restart-file over data-file(s)
            lmp.command("read_restart {}".format(args.restart))
        elif args.fdat:

            for cur_data_nr, cur_lmpdat in enumerate(args.fdat):
                if cur_data_nr < 1:
                    lmp.command("read_data {}".format(cur_lmpdat))
                else:
                    # append atoms from other data-files
                    lmp.command("read_data {} add append".format(cur_lmpdat))

        else:
            raise RuntimeError("Neither restart nor data file defined!")

    # not first run, atom number changed
    elif new_run:
        lmp.command("read_data restart.lmpdat")
    else:
        lmp.command("read_restart resume.lmprst")

    # read file with long-range interactions (pair-coefficients)
    lmp.file(args.fpc)

    # ==========================================================================#
    # INTEGRATION =============================================================#
    # ==========================================================================#
    lmp.command("fix nvt_1 all nvt temp {} {} {}".format(Tstart, Tstop, 100))

    # ==========================================================================#
    # OUTPUT-SETTINGS =========================================================#
    # ==========================================================================#
    # * Restart-/dcd-filenames
    f_data = "{0}-{1}.lmpdat".format(args.o, cur_iteration + 1)
    f_restart = "{0}-{1}.lmprst".format(args.o, cur_iteration)
    f_dcd = "{0}-{1}.dcd".format(args.o, cur_iteration)

    # * Thermodynamic-logging
    lmp.command("thermo {}".format(thermsteps))
    Thermargs = [
        "step",
        "temp",
        "press",
        "vol",
        "density",
        "cella",
        "cellb",
        "cellc",
        "cellalpha",
        "cellbeta",
        "cellgamma",
        "atoms",
        "bonds",
        "angles",
        "etotal",
        "pe",
        "evdwl",
        "ecoul",
        "ebond",
        "eangle",
        "edihed",
        "eimp",
        "enthalpy",
    ]
    lmp.command("thermo_style custom " + " ".join(Thermargs))
    # modify output, thermodynamics, etc.
    lmp.command("thermo_modify lost warn flush yes")

    # ==========================================================================#
    # SIMULATE ================================================================#
    # ==========================================================================#
    # Set or change the velocities of a group of atoms (here: all)
    lmp.command("velocity all create {} 483806 rot yes dist gaussian".format(Tstop))

    # write trajectory every Dsteps
    lmp.command("dump ID_DCD1 all dcd {} {}".format(Dsteps, f_dcd))

    # ==========================================================================#
    # EQUILIBRATE =============================================================#
    # ==========================================================================#
    lmp.command("run {}".format(args.steps))

    # ==========================================================================#
    # CHECK COORDINATES AND DELETE MOLECULES BEYOND BORDER ====================#
    # ==========================================================================#
    molecule_ids = lmp.gather_atoms("molecule", 0, 1)
    lmp_coords = lmp.gather_atoms("x", 1, 3)  # list of coordinates as floats

    # Serial part; find molecules with largest z-coordinates
    if rank == 0:
        ids_to_delete = []
        coords = list(lmp_coords)

        for z_coord in coords[2::3]:
            if z_coord > zmax1 or z_coord < zmax2:
                del_molecule_id = coords.index(z_coord) / 3
                ids_to_delete.append(molecule_ids[del_molecule_id])

        # remove molecule-id-duplicates
        ids_to_delete = set(ids_to_delete)
    else:
        # in the meantime other nodes just store a placeholder
        ids_to_delete = None

    # Back to parallel, broadcast outcome from root-node to all remaining nodes
    ids_to_delete = comm.bcast(ids_to_delete, root=0)

    # only do a full restart when there are more than XY molecules to delete
    if len(ids_to_delete) > 10:
        for mol_id in ids_to_delete:
            # create group with current mol-id
            lmp.command("group remove_molecule molecule {}".format(mol_id))
            # delete all atoms from group
            lmp.command("delete_atoms group remove_molecule mol yes")
            # delete group
            lmp.command("group remove_molecule delete")

        # ======================================================================#
        # WRITE NEW DATA (TOPOLOGY)-, RESTART- AND DCD-FILE ===================#
        # (ATOM-NR's do not match anymore)                  ===================#
        # ======================================================================#
        # get new topology
        # for restarting the simulation, gets overwritten each iteration
        lmp.command("write_data restart.lmpdat")
        # *********** write data for the moment, later on, convert restart to psf
        lmp.command(
            "write_data {}-{}.lmpdat".format(args.o, cur_iteration)
        )  # for the dcd

        if rank == 0:
            # data for next iteration
            lmga.refix_lmp_data("restart.lmpdat")

        lmp.command("undump ID_DCD1")
        new_run = True
    else:
        lmp.command("write_restart resume.lmprst")
        new_run = False

    # * cleaning
    if rank == 0:
        sl.move("log.lammps", "{0}-{1}.lmplog".format(args.o, cur_iteration))

    del (molecule_ids, lmp_coords, ids_to_delete)
    lmp.close()

    # first
    first_run = False
