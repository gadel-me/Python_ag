from __future__ import print_function, division
import os
import copy
import shutil as sl
import re
import argparse
import math
import pdb
import time
import numpy as np
from natsort import natsorted
#import itertools as it
import scipy.stats
from mpi4py import MPI
from lammps import lammps, PyLammps
import Transformations as cgt
import md_elements as mde
import md_box as mdb
import ag_unify_md as agum
import ag_geometry as agm
import ag_lammps as aglmp
import ag_unify_log as agul
import ag_vectalg as agv

#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator

#==============================================================================#
# Helper functions and variables
#==============================================================================#
initial_velocity = "velocity all create {} 8324798 rot yes dist gaussian"
_sqrt_3 = math.sqrt(3)


def lmp_thermostat(lmp, logsteps):
    """
    Log thermodynamic data.
    """
    thermargs = ["step", "temp", "press", "vol", "density",
                 "cella", "cellb", "cellc", "cellalpha", "cellbeta", "cellgamma",
                 "etotal", "pe", "evdwl", "ecoul", "ebond", "eangle", "edihed",
                 "eimp", "enthalpy"]
    lmp.command("thermo_style custom " + " ".join(thermargs))
    lmp.command("thermo_modify lost warn flush yes")
    #lmp.command("thermo_modify line multi format float %g")
    lmp.command("thermo {}".format(logsteps))


def mkdir(folder):
    """
    Create folder or skip creation if it already exists.

    Create a folder, if it exists, do nothing.
    """
    try:
        os.mkdir(folder, 0755)
    except OSError:
        print("***Info: Folder {} already exists!".format(folder))


def _read_system(lmpdat, dcd=None, frame_idx=-1):
    """
    Read a lammps data file and a dcd file on top (optional).

    This function simplifies the process of reading a lammps data file and loading
    a dcd file on top of it for further coordinates. It returns a Universe object
    which has many methods to manipulate it with.

    Parameters
    ----------
    lmpdat : str
        Name of the lammps data file

    dcd : str, optional
        Name of the dcd file with the same amounts of atoms as the lammps data
        file

    frame_idx : int, optional
        Index of the frame to use from the dcd file.

    Returns
    -------
    md_sys : LmpStuff object
        An object of LmpStuff which can be further processed.

    """
    # read output from quenching
    md_sys = aglmp.LmpStuff()
    md_sys.read_lmpdat(lmpdat)

    if dcd is not None:
        md_sys.import_dcd(dcd)
        # since we are only interested in one frame, delete all others
        md_sys.ts_coords = []
        md_sys.ts_boxes = []

        # enable reading the last frame with negative indexing
        if frame_idx == -1:
            md_sys.read_frames(frame=frame_idx - 1, to_frame=frame_idx)
        else:
            md_sys.read_frames(frame=frame_idx, to_frame=frame_idx + 1)

        md_sys.close_dcd()

    return md_sys


def check_aggregate(mdsys, frame_id=-1, atm_atm_dist=4, unwrap=False, debug=False):
    """
    Check if several molecules form an aggregate.

    #TODO Buggy, does not recognize aggregates properly when they form a chain?

    At the moment this works only for orthogonal cells (reason: sub-cell length)
    Check if the aggregate did not get dissolved in the process. This is achieved
    by calculating the center of geometry for each molecule and subsequent com-
    parison of the distance between each center of geometry. The smallest distance
    must not be larger than four times the radius of the biggest molecule in the
    system.

    Input:
        > mdsys             ag_unify_md.Unification; system output from 'SYS-PREPARATION'-step
        > frame_id          int; frame to process
        > atm_atm_dist      float; radius around an atom in which another atom
                            should be positioned
        > debug             boolean; True if further output should be given,
                            default=False

    Output:
        > aggregate_ok      boolean; True, if aggregate is still intact, False if
                            a part of it drifted away
    """
    aggregate_ok = False

    # lammps may internally have the coordinates wrapped -> unwrap cell when
    # necessary (e.g. coordinates were taken from restart)
    if unwrap is True:
        print("***Check Aggregate Info: Unwrapping cell")
        mdsys.unwrap_cell(frame_id)

    # maximal possible distance between two atoms
    cube_side = atm_atm_dist / _sqrt_3  # -> atm_atm_dist/(2*math.sqrt(3)) * 2

    if debug is True:
        print("***Check Aggregate Info: Cube side {}".format(cube_side))

    mdsys.create_linked_cells(frame_id, rcut_a=cube_side, rcut_b=cube_side,
                              rcut_c=cube_side)

    # include same molecule or it gets missing in our aggregates set
    close_atoms, aggregates = mdsys.chk_atm_dist(frame_id,
                                                 min_dist=atm_atm_dist,
                                                 exclude_same_molecule=False,
                                                 get_aggregates=True,
                                                 debug=False)

    if debug is True:
        print("***Check Aggregate Info: Number of cubes in direction a {}".format(mdsys.ts_lnk_cls[frame_id].ra))
        print("***Check Aggregate Info: Number of cubes in direction b {}".format(mdsys.ts_lnk_cls[frame_id].rb))
        print("***Check Aggregate Info: Number of cubes in direction c {}".format(mdsys.ts_lnk_cls[frame_id].rc))
        print("***Group IDs of all aggregates: {}".format(aggregates))
        print("***IDs of close atoms: {}".format(" ".join([str(i) for i in close_atoms])))

    # aggregate is only o.k. if all molecules are part of it
    if len(aggregates) == 1 and len(aggregates[0]) == len(mdsys.molecules):
        aggregate_ok = True

        if debug is True:
            print("***Check Aggregate Info: Aggregate looks fine!")

    return aggregate_ok


# System Preparation ==========================================================#

def _rot_matrix():
    """
    Create an arbitrary rotational matrix.

    Returns
    -------
    m_arb_xyz : array of arrays
        The arbitrary rotational matrix

    """
    # rotation matrix, consisting of 3 arbitrary rotations A, B and C
    m_arb_rot_x = agm.arb_rot_matrix([1, 0, 0])  # A
    m_arb_rot_y = agm.arb_rot_matrix([0, 1, 0])  # B
    m_arb_rot_z = agm.arb_rot_matrix([0, 0, 1])  # C
    # multiply all three matrices; ABC=A(BC)=(AB)C
    m_arb_xy = np.matmul(m_arb_rot_x, m_arb_rot_y)  # AB
    m_arb_xyz = np.matmul(m_arb_xy, m_arb_rot_z)  # ABC
    return m_arb_xyz


def _rotate_sys(md_sys):
    """
    Rotate system arbitrarily.
    """
    m_arb_xyz = _rot_matrix()
    sys_add_all_atms = range(len(md_sys.atoms))
    md_sys.mm_atm_coords(-1, m_arb_xyz, False, *sys_add_all_atms)


def _shift_sys(md_sys, radius, radius_buffer=1):
    """
    """
    # enlarge radius so atoms do not clash
    radius += radius_buffer
    rn_pos = agm.points_on_sphere(npoints=1, ndim=3, radius=radius)[0]
    mx_trans = cgt.translation_matrix(rn_pos)
    sys_add_all_atms = range(len(md_sys.atoms))
    md_sys.mm_atm_coords(-1, mx_trans, False, *sys_add_all_atms)


def _create_new_box(md_sys):
    """
    """
    # add new orthogonal box
    md_sys.ts_boxes = []  # delete all previous unneeded boxes
    box_diameter = md_sys.get_system_radius(-1) + 30
    pi_2 = math.pi / 2
    new_box = mdb.Box(boxtype="lattice", ltc_a=box_diameter, ltc_b=box_diameter,
                      ltc_c=box_diameter, ltc_alpha=pi_2, ltc_beta=pi_2,
                      ltc_gamma=pi_2)
    md_sys.ts_boxes.append(new_box)
    md_sys.ts_boxes[0].box_lat2lmp(triclinic=False)


def _check_sys(md_sys, atm_idx_max):
    """
    """
    # check interatomic distances
    md_sys.create_linked_cells(-1)
    close_atms = md_sys.chk_atm_dist(frame_id=-1, min_dist=1.2)
    close_atms = sorted(close_atms)

    # check if new atoms were placed too near to the main system by comparing
    # the atom indices (atom indices lower than those from the original main
    # system are not counted)
    try:
        return close_atms[-1] > atm_idx_max
    except IndexError:
        return True


def sysprep(lmpdat_out, lmpdat_main, lmpdat_add, dcd_add=None,
            frame_idx=-1):
    """
    Prepare the system for the next docking step.

    Since the kawska-zahn approach adds an agglomerate of molecules each stage,
    the system must be prepared for lammps beforehand. Therefor a sphere with
    radius r1 is created around the main system as well as a sphere with radius
    r2 around the agglomerate to add. Both systems are combined, checked for clashes
    and eventually a new lammps data file is written which can further be
    utilized. If

    Parameters
    ----------
    lmpdat_out : str
        Name of merged lammps data file to write

    lmpdat_main : str
        Name of the lammps data file to add another system to

    lmpdat_add : str
        Name of the lammps data file which is added to the main molecular
        system

    Returns
    -------
    success : bool
        True if successful, False otherwise.

    Writes a new lammps data file called sysprep_out_index

    """
    # read and transpose the main sys to the origin
    main_sys = _read_system(lmpdat_main)
    main_sys.transpose_by_cog(-1, [0, 0, 0], copy=False)
    _natoms = len(main_sys.atoms)

    # read and transpose the add sys to the origin
    add_sys = _read_system(lmpdat_add, dcd_add, frame_idx)
    add_sys.transpose_by_cog(-1, [0, 0, 0], copy=False)

    # rotate add sys
    _rotate_sys(add_sys)

    # shift add sys to sphere around main sys
    main_sys_radius = main_sys.get_system_radius(-1)
    add_sys_radius = add_sys.get_system_radius(-1)
    kwz_radius = main_sys_radius + add_sys_radius

    _shift_sys(add_sys, kwz_radius)

    # merge both systems
    main_sys.extend_universe(add_sys, mode="merge")
    _create_new_box(main_sys)

    # group atoms by bonds
    main_sys.fetch_molecules_by_bonds()
    main_sys.mols_to_grps()

    success = _check_sys(main_sys, _natoms)

    # write an output lammps data only if everything worked out
    if success is True:
        # write new data file
        main_sys.change_indices(incr=1, mode="increase")
        main_sys.write_lmpdat(lmpdat_out, frame_id=0, title="System ready" +
                              "for docking", cgcmm=True)
    return success


# Quenching ===================================================================#
def quench(lmpdat_files, settings_files, temps, steps, quench_output, gpu):
    """
    """
    main_sys = _read_system(lmpdat_files[0])
    natoms_main_sys = len(main_sys.atoms)
    del main_sys

    quench_lmp = lammps()
    quench_pylmp = PyLammps(ptr=quench_lmp)
    quench_lmp.command("log {} append".format(quench_output[0]))

    # load gpu package
    if gpu:
        # neighbor list building on gpu (currently) leads to segmentation faults
        #quench_lmp.command("package gpu 1 neigh yes")
        quench_lmp.command("package gpu 1 neigh no")  # DEBUGGING
        quench_lmp.command("suffix gpu")

    # load general settings (units, boundary, dimension, atom_style, etc. pp.)
    quench_lmp.file(settings_files[0])

    # change box type to not be periodic
    quench_lmp.command("boundary f f f")

    # read data and create an initial velocity or read the restart file from
    # a previous run
    if os.path.isfile(quench_output[2]) is True:
        quench_lmp.command("read_restart {}".format(quench_output[2]))
    else:
        quench_lmp.command("read_data {}".format(lmpdat_files[1]))
        #quench_lmp.command(initial_velocity.format(temps[0]))

    # read pair coefficients from file if provided
    if settings_files[1] is not None:
        quench_lmp.file(settings_files[1])

    lmp_thermostat(quench_lmp, steps[0])
    # define the atoms that may move during the simulation
    quench_lmp.command("group grp_add_sys id > {}".format(natoms_main_sys))
    quench_lmp.command("group grp_main_sys id <= {}".format(natoms_main_sys))
    quench_lmp.command("fix freeze grp_main_sys setforce {0} {0} {0}".format(0.0))

    # barostatting, thermostatting
    #quench_lmp.command("fix thermostat_barostat grp_add_sys nvt temp {0} {1} 0.1".format(*temps))

    # set an additional push to the added atoms that should be docked
    # towards the origin at (0/0/0)
    # (prevents losing atoms due to being localized outside the cutoff)
    if rank == 0:
        prep_sys = _read_system(lmpdat_files[1])
        cog = agm.get_cog(prep_sys.ts_coords[-1][natoms_main_sys + 1:])
        cog /= np.linalg.norm(cog, axis=0)  # unit vector
        cog *= -50
    else:
        cog = None

    cog = comm.bcast(cog, 0)
    quench_lmp.command(("fix force grp_add_sys addforce {0} {1} {2} every 1000").format(*cog))

    # trajectory
    quench_lmp.command("dump trajectory all dcd {} {}".format(steps[0], quench_output[1]))
    quench_lmp.command("dump_modify trajectory unwrap yes")

    # intermediate restart files
    #tmp1_rst  = quench_dir + "tmp-1-quench_out_{}.lmprst".format(curcycle)
    quench_lmp.command("restart {0} {1} {1}".format(steps[0] * 10, quench_output[2]))

    # 1st run with a little push towards the aggregate
    #if os.path.isfile(quench_output[2]) is False:
    #    quench_lmp.command("run {}".format(20000))

    # remove artificial pushing force and minimize the docked complex
    #quench_lmp.command("unfix group_loose_addforce")
    #quench_lmp.command("run {}".format(50000))
    quench_lmp.command("min_style fire")
    quench_lmp.command("minimize 1.0e-5 1.0e-8 1000 1000000")
    quench_lmp.command("min_style cg")
    quench_lmp.command("minimize 1.0e-5 1.0e-8 1000 1000000")

    # write restart file
    quench_lmp.command("undump trajectory")
    quench_lmp.command("unfix freeze")
    quench_lmp.command("unfix force")
    #quench_lmp.command("unfix thermostat_barostat")
    quench_lmp.command("reset_timestep 0")
    quench_lmp.command("write_restart {}".format(quench_output[2]))
