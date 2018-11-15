#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import numpy as np
import argparse
import os
import re
from collections import OrderedDict
from lammps import lammps
from mpi4py import MPI
import ag_unify_md as agum
import ag_unify_log as agul
import ag_geometry as agg
import pdb
import time


#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator

#==============================================================================#
# Helping functions
#==============================================================================#
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def get_scanned_geometry(gau_log):
    """
    Bla.

    Search for keyword scan in gaussian output file and return type of geometry
    (e.g. bond, angle, dihedral) that was scanned.
    """
    with open(gau_log) as opened_gau_log:
        line = opened_gau_log.readline()
        while line != "":
            if "Initial Parameters" in line:
                while "Scan" not in line:
                    line = opened_gau_log.readline()
                else:
                    split_line = line.split()
                    scanned_geometry = split_line[2]
                    geometry_value = float(split_line[3])
                    break

            line = opened_gau_log.readline()

    return (scanned_geometry, geometry_value)


def get_geometry_by_key(key_string):
    """
    Geometry type by string length.

    Get type of the geometry (bond/angle/dihedral) by the number of atom indices.
    """
    split_key = key_string.split()

    if len(split_key) == 2:
        cur_geometry = "bond"
    elif len(split_key) == 3:
        cur_geometry = "angle"
    elif len(split_key) == 4:
        cur_geometry = "dihedral"
    else:
        raise IOError("Number of atom indices is wrong")

    return cur_geometry


def built_restrain_string(indices_and_values, k_start=0.0, k_stop=5000):
    """
    Bla.

    Define a string for fix restrain by a dictionary with given atom indices
    as strings and their according lengths/angles as values.
    """
    restrain_string = "fix REST all restrain "

    for key, value in indices_and_values.iteritems():
        cur_geometry = get_geometry_by_key(key)
        restrain_string += "{} {} {} {} {}".format(cur_geometry, key, k_start, k_stop, value)

    return restrain_string


def scan(lmpdat, output, indices_and_values, temp=(600, 0), k=(0.0, 200.0)):
    """
    Scan a bond, angle or dihedral.

    Carry out a force field calculation for an atom, a bond or a dihedral using
    the fix restrain command.

    lmpdat              str; lammps data file
    indices_and_values  dict; str with geometry indices is the key and value
                        is the value of the geometry
                        e.g. {"3 4 5 12": 180.05, "4 5 2 6 1": 168.90, ...}
    k_start             float; starting value for force constant k
    k_stop              float; stopping value for force constant k
    output              str; appendix to output files

    Sources:    (https://lammps.sandia.gov/threads/msg17271.html)
                https://lammps.sandia.gov/doc/fix_restrain.html
    """
    save_step   = 50000
    anneal_step = 750000
    quench_step = 500000
    #anneal_step = 100000
    #quench_step = 100000
    thermargs   = ["step", "temp", "pe", "eangle", "edihed", "eimp", "evdwl", "ecoul", "ebond", "enthalpy"]

    # split world communicator into n partitions and run lammps only on that
    # specific one
    split = comm.Split(rank, key=0)
    lmp = lammps(comm=split)

    lmp.command("log {}.lmplog".format(output))
    lmp.command("units metal")
    lmp.command("boundary p p p")
    lmp.command("dimension 3")
    lmp.command("atom_style full")
    lmp.command("pair_style lj/cut/coul/cut 30.0")
    lmp.command("bond_style harmonic")
    lmp.command("angle_style harmonic")
    lmp.command("dihedral_style charmm")
    lmp.command("improper_style cvff")
    lmp.command("special_bonds amber")
    lmp.command("timestep 0.001")

    # data and thermo
    lmp.command("read_data {}".format(lmpdat))

    # make lammps calculate the value of the entity (bond, angle, dihedral)
    lmp.command("thermo_style custom " + " ".join(thermargs))
    lmp.command("thermo_modify lost warn flush yes")
    lmp.command("thermo {}".format(save_step))

    lmp.command("fix MOMENT all momentum 100 linear 1 1 1 angular")
    lmp.command("dump DUMP all dcd {} {}.dcd".format(save_step, output))

    # annealing -> achieve wanted angle (may be omitted?)
    print(bcolors.OKGREEN + "Annealing\n" + bcolors.ENDC)

    ###########################################################################
    # TESTING
    #indices_and_values = {"25 8 9 10" : 90.0}
    #indices_and_values = {"11 7 8 25" : -4.0}
    #indices_and_values = {"25 8 7 11" : 3.1415}
    #lmp.command("group entity_atoms id {}".format(indices_and_values.keys()[0]))
    #lmp.command("compute 1 entity_atoms property/local dtype datom1 datom2 datom3 datom4")
    #lmp.command("compute 2 entity_atoms dihedral/local phi")
    #lmp.command("dump 1 all local 1000 {}.dump index c_2".format(output))
    ###########################################################################

    restrain_anneal = built_restrain_string(indices_and_values, k[0], k[1])
    lmp.command("velocity all create {} 8675309 mom yes rot yes dist gaussian".format(temp[0]))
    lmp.command("fix NVE all nve")
    lmp.command("fix TFIX all langevin {} {} 100 24601".format(temp[0], temp[1]))
    #lmp.command("fix TFIX all temp/berendsen {} {} 5".format(temp[0], temp[1]))
    #lmp.command("fix TFIX all nvt temp {} {} 0.1".format(temp[0], temp[1]))
    lmp.command(restrain_anneal)
    lmp.command("fix_modify REST energy yes")
    lmp.command("run {}".format(anneal_step))

    # quenching
    print(bcolors.OKGREEN + "Quenching\n" + bcolors.ENDC)
    restrain_quench = built_restrain_string(indices_and_values, k[1], k[1])
    lmp.command("fix TFIX all langevin {} {} 100 24601".format(temp[1], temp[1]))
    #lmp.command("fix TFIX all langevin {} {} 10 24601".format(temp[1], 0))
    #lmp.command("fix TFIX all nvt temp {} {} 0.1".format(temp[1], temp[1]))
    lmp.command(restrain_quench)
    lmp.command("fix_modify REST energy yes")
    lmp.command("run {}".format(quench_step))

    # sanity check for convergence
    print(bcolors.OKGREEN + "Minimization\n" + bcolors.ENDC)
    # output
    lmp.command("minimize 1e-6 1e-9 2000000 100000")
    # report unrestrained energies, single point energy
    lmp.command("unfix REST")
    lmp.command("run 0")
    lmp.close()


def measure_geometry(lmpdat, dcd, atm_ids):
    """
    Bla.

    Measure the actual bond, angle, dihedral achieved by restraining it.
    """
    sys_md = agum.Unification()
    sys_md.read_lmpdat(lmpdat)
    sys_md.import_dcd(dcd)
    sys_md.read_frames(frame=-2,
                       to_frame=-1,
                       frame_by="index")  # read only the last frame
    # decrement atom indices since for sys_md they start with 0, whereas
    # for lammps and gaussian they start with 1
    atm_idxs = [int(i) - 1 for i in atm_ids.split(" ")]

    # dummy value, good for debugging
    geometry_value = None

    # bond
    if len(atm_idxs) == 2:
        geometry_value = np.linalg.norm(sys_md.ts_coords[-1][atm_idxs[0]] -
                                        sys_md.ts_coords[-1][atm_idxs[1]])
    # angle
    elif len(atm_idxs) == 3:
        geometry_value = agg.get_angle(sys_md.ts_coords[-1][atm_idxs[0]],
                                       sys_md.ts_coords[-1][atm_idxs[1]],
                                       sys_md.ts_coords[-1][atm_idxs[2]])
    # dihedral
    elif len(atm_idxs) == 4:
        geometry_value = agg.new_dihedral([sys_md.ts_coords[-1][atm_idxs[0]],
                                           sys_md.ts_coords[-1][atm_idxs[1]],
                                           sys_md.ts_coords[-1][atm_idxs[2]],
                                           sys_md.ts_coords[-1][atm_idxs[3]]])
    else:
        raise IOError("Wrong number of atom ids!")

    if len(atm_idxs) > 2:
        geometry_value = np.degrees(geometry_value)

    return geometry_value


def get_entities(gau_log, geom_entity):
    """
    Bla.

    Extract all entities from gaussian output file.
    """
    # pull all other entities (bonds, angles, dihedrals) from gaussian output
    sys_gau = agum.Unification()
    sys_gau.read_gau_log(gau_log)
    entities = sys_gau.gaussian_other_info[geom_entity]
    return entities


def get_last_pe_value(lmplog):
    """
    Bla.

    Pull energy value from lammps calculation.
    """
    log_data = agul.LogUnification()
    log_data.read_lmplog(lmplog)
    return log_data.data[-1]["PotEng"][-1]


def get_entity(dictionary):
    """
    """
    key_length = len(dictionary.keys()[0])
    if key_length == 2:
        return "Bond [Angstrom]"
    elif key_length == 3 or key_length == 4:
        return "Angle [Degrees]"
    else:
        return "Unknown"


def write_energies_file_header(filename):
    """
    Bla.

    Write header for energies file.
    """
    with open(filename, "w") as opened_filename:
        opened_filename.write("{:>20s} {:>20s}\n".format("Angle [Degrees]",
                                                         "Energy [eV]"))


def md_from_ab_initio(gau_log, lmpdat, temp=(600, 0), k=[0.0, None],
                      energy_file_out="defaultname", output_idx=0):
    """
    Calculate md energy.

    Extract all necessary info from a gaussian log file, perform a simulated
    annealing simulation to get a relaxed structure and perform a minimization
    afterwards. This routine follows the path of the lowest energy, keeping
    the bond, angle or dihedral of interest the closest to that outcome of the
    ab initio calculation. Increasing the temperature or in-/decreasing the
    force constant k may be necessary.

    See: https://lammps.sandia.gov/doc/fix_restrain.html

    gau_log     str; gaussian log file
    lmpdat      str; lammps data file
    output      str; output-prefix for the data being written
    temp        tuple; starting and finishing temperature for the annealing
                procedure
    k           tuple; starting and stopping force constants for the geometric
                entity being examined
    """

    # grab scanned entity from ab initio output and involved atom ids
    geom_entity, _ = get_scanned_geometry(gau_log)
    ids_geom_enitity = re.findall(r"\d+", geom_entity)

    # split single mds over the cores available
    entities = get_entities(gau_log, geom_entity)

    if rank == 0:
        if not os.path.isfile(energy_file_out):
            write_energies_file_header(energy_file_out)

    if geom_entity.startswith("R"):
        if k[1] is None:
            k[1] = 1200
    elif geom_entity.startswith("A"):
        if k[1] is None:
            k[1] = 300
    # shift phase by 180 degrees as defined in lamps manual
    # See: https://lammps.sandia.gov/doc/fix_restrain.html, "dihedral"
    elif geom_entity.startswith("D"):
        if k[1] is None:
            k[1] = 80
    else:
        raise Warning("Entry seems odd")

    for task, cur_geom_value in enumerate(entities):

        # parallelization; each rank does this loop and skips it, if it is not
        # its turn
        if task % size != rank:
            continue

        # shift phase by 180 degrees as defined in lammps manual
        # See: https://lammps.sandia.gov/doc/fix_restrain.html, "dihedral"
        if geom_entity.startswith("D"):
            # shift dihedrals by 180 degrees (dont know why)
            cur_geom_value += 180

        # define a dictionary with atom ids as key and the current
        # geometry value as value
        cur_geom_atm_ids = " ".join(ids_geom_enitity)
        cur_geom_atm_ids_geom_value = {cur_geom_atm_ids: cur_geom_value}

        # define appendix for all files
        output_appendix = "{}_{}_{}".format(cur_geom_atm_ids.replace(" ", "_"),
                                            task, output_idx)

        # skip calculations that were already carried out
        if os.path.isfile(output_appendix + ".lmplog"):
            print(output_appendix + ".lmplog" + " already exists!")
            continue

        # annealing with quenching and minimization using lammps
        scan(lmpdat, output_appendix, cur_geom_atm_ids_geom_value, temp, k)

        # since minimized md-structure != minimized ab initio structure,
        # the real value of the geometry of interest must be derived from the
        # last frame of the md simulation
        current_geom_val = measure_geometry(lmpdat, output_appendix + ".dcd",
                                            cur_geom_atm_ids)

        # pull energy from md run
        last_pe_value = get_last_pe_value(output_appendix + ".lmplog")

        # pre defined string for each row in results file
        resume_file_row = "{:> 20.4f} {:> 20.8f}\n"

        with open(energy_file_out, "a") as opened_resume_file:
            opened_resume_file.write(resume_file_row.format(current_geom_val,
                                                            last_pe_value))

        print("{} finished".format(rank))
        #time.sleep(5)


def norm_energy(energy_file_in, energy_file_out):
    """
    Bla.

    Norm energies by smallest value and sort by angle.
    """
    keys = []
    original_values = []

    with open(energy_file_in) as opened_energy_file:
        # skip header
        opened_energy_file.readline()
        line = opened_energy_file.readline()

        while line != "":
            split_line = line.split()
            keys.append(float(split_line[0]))
            original_values.append(float(split_line[1]))
            line = opened_energy_file.readline()

    min_value = min(original_values)
    normed_values = []

    for val in original_values:
        val -= min_value
        normed_values.append(val)

    keys_and_values = dict(zip(keys, normed_values))
    keys_and_values = OrderedDict(sorted(keys_and_values.items()))

    if not os.path.isfile(energy_file_out):
        write_energies_file_header(energy_file_out)

    #with open(energy_file_out, "a") as opened_energy_file:
    with open(energy_file_out, "a") as opened_energy_file:

        for key, value in keys_and_values.iteritems():
            opened_energy_file.write("{:> 20.8f} {:> 20.8f}\n".format(key, value))


#==============================================================================#
# User Input
#==============================================================================#
if __name__ == "__main__":
    if rank == 0:
        parser = argparse.ArgumentParser()
        parser.add_argument("lmpdat",
                            help="Lammps' data-file")

        parser.add_argument("-gau_logs",
                            nargs="*",
                            help="Gaussian log-file")

        parser.add_argument("-geometries",
                            nargs="*",
                            help="Atom indices, where each string is one set of one geometry")

        parser.add_argument("-out",
                            default="DEFAULTNAME",
                            help="Name of energy files.")

        args = parser.parse_args()
    else:
        args = None

    args = comm.bcast(args, 0)

    #==============================================================================#
    # do restrained md for geometry scanned in gaussian
    #==============================================================================#
    output_file = "{}_md.txt".format(args.out)

    if not os.path.isfile(output_file):
        for gau_file_idx, cur_gau_log in enumerate(args.gau_logs):
            # Use k=80 for dihedrals and k=200 for angles and k=1200 bonds
            #md_from_ab_initio(cur_gau_log, args.lmpdat, energy_file_out=output_file, output_idx=gau_file_idx, temp=(600, 0), k=(0.0, 1200.0))
            md_from_ab_initio(cur_gau_log, args.lmpdat, energy_file_out=output_file, output_idx=gau_file_idx)

    # wait for all ranks to finish
    time.sleep(2)
    print("{} is done".format(rank))

    # norm energies
    if rank == 0:
        norm_energy(output_file, "{}_md_normed.txt".format(args.out))
