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
import ag_lammps as aglmp
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

# split world communicator into n partitions and run lammps only on that
# specific one
split = comm.Split(rank, key=0)

#==============================================================================#
# Helping functions
#==============================================================================#
class bcolors:
    header = '\033[95m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    endc = '\033[0m'
    bold = '\033[1m'
    underline = '\033[4m'


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


def compile_restrain_string(indices_and_values, force_constants, hold=0):
    """
    Compile the lammps string to restrain certain geometry entities.

    Parameters
    ----------
    indices_and_values : dict
        atom ids and according entity values, e.g. {"1 2 3 4": 120, ...}

    force_constants : list of lists
        starting and ending values of the force constants to use

    hold : 0 or 1
        assign the force constant to ending and start (1) or just
        for the ending (0) of the restrain string

    Define a string for fix restrain by a dictionary with given atom indices
    as strings and their according lengths/angles as values.
    """
    # short check if function is used the right way
    assert hold in (0, 1)
    #pdb.set_trace()

    restrain_string = "fix REST all restrain "

    for index, (key, value) in enumerate(indices_and_values.iteritems()):
        k_start = force_constants[index][hold]
        k_stop = force_constants[index][1]
        cur_geometry = get_geometry_by_key(key)
        restrain_string += "{} {} {} {} {} ".format(cur_geometry, key, k_start, k_stop, value)

    return restrain_string


def null_entity(lmp, lmpdat, indices_and_values, key_index=0):
    """
    Does not work atm. Instead, read and write a new data file with dummy
    bond-/angle- or dihedraĺ-coeffs.

    Compile a bond-/angle-/dihedral-coeff string for lammps and alter the lammps instance.

    The 'coeff-string' has only values of 0 which means it gets omitted.
    Needed to switch off the energy contribution of a given dihedral in the
    system.

    Parameters
    ----------
    lmp : lammps.lammps instance
        the instance of lammps which will be altered
    lmpdat : str
        lammps data file that is to be read
    indices_and_values : dict
        atom ids and according entity values, e.g. {"1 2 3 4": 120, ...}
    key_index : int
        index of indices_and_values to be processed
    """
    #lmp_sys = aglmp.read_lmpdat(lmpdat)
    #num_geom_types = None
    #key = indices_and_values.keys()[key_index]
    #cur_geometry = get_geometry_by_key(key)
    #null_coeff = None
#
    #if cur_geometry == "bond":
    #    num_geom_types = len(lmp_sys.bnd_types)
    #    null_coeff_id = num_geom_types + 1
    #    null_coeff = "bond_coeff {} 0.0 0.0".format(null_coeff_id)
    #elif cur_geometry == "angle":
    #    num_geom_types = len(lmp_sys.ang_types)
    #    null_coeff_id = num_geom_types + 1
    #    null_coeff = "angle_coeff {} 0.0 0.0".format(null_coeff_id)
    #elif cur_geometry == "dihedral":
    #    num_geom_types = len(lmp_sys.dih_types)
    #    null_coeff_id = num_geom_types + 1
    #    null_coeff = "dihedral_coeff {} 0.0 1 0 0.5".format(null_coeff_id-1)
    #else:
    #    raise Warning("***Warning: Something went wrong!")
#
    #pdb.set_trace()
    #lmp.command(null_coeff)
    #lmp.command("set atom {} {} {}".format(key, cur_geometry, null_coeff_id))


def scan(lmpdat, output, indices_and_values, ks, temps=(600, 0), omit_entity=False):
    """
    Scan a bond, angle or dihedral.

    Carry out a force field calculation for an atom, a bond or a dihedral using
    the fix restrain command.

    Parameters
    ----------
    lmpdat : str
        lammps data file
    indices_and_values : dict
        geometry indices is the key and value is the value of the geometry
        e.g. {"3 4 5 12": 180.05, "4 5 2 6 1": 168.90, ...}
    output : str
        appendix to output files

    Sources:    (https://lammps.sandia.gov/threads/msg17271.html)
                https://lammps.sandia.gov/doc/fix_restrain.html

    """
    save_step = 50000
    anneal_step = 750000
    quench_step = 500000
    thermargs = ["step", "temp", "pe", "eangle", "edihed", "eimp", "evdwl",
                 "ecoul", "ebond", "enthalpy"]

    # split world communicator into n partitions and run lammps only on that specific one
    lmp = lammps(comm=split, cmdargs=["-screen", "lmp_out.txt"])

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
    print(bcolors.red + "Annealing on rank {}: ".format(rank) + output + bcolors.endc)
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

    lmp.command("velocity all create {} 8675309 mom yes rot yes dist gaussian".format(temps[0]))
    lmp.command("fix NVE all nve")
    lmp.command("fix TFIX all langevin {} {} 100 24601".format(temps[0], temps[1]))
    restrain_anneal = compile_restrain_string(indices_and_values, ks, hold=0)
    lmp.command(restrain_anneal)
    lmp.command("fix_modify REST energy yes")

    try:
        lmp.command("run {}".format(anneal_step))
    except Exception:
        print("***Error: Simulation crashed (annealing)! Force constants too high?")
        MPI.COMM_WORLD.Abort()

    # quenching
    print(bcolors.yellow + "Quenching on rank {}: ".format(rank) + output + bcolors.endc)
    restrain_quench = compile_restrain_string(indices_and_values, ks, hold=1)
    lmp.command("fix TFIX all langevin {} {} 100 24601".format(temps[1], temps[1]))
    lmp.command(restrain_quench)
    lmp.command("fix_modify REST energy yes")

    try:
        lmp.command("run {}".format(quench_step))
    except:
        print("***Error: Simulation crashed (quenching)! Force constants too high?")
        MPI.COMM_WORLD.Abort()

    # sanity check for convergence
    print(bcolors.green + "Minimization on rank {}: ".format(rank) + output + bcolors.endc)

    try:
        lmp.command("minimize 1e-6 1e-9 2000000 100000")
    except:
        print("***Error:  Simulation crashed (minimization)!")
        MPI.COMM_WORLD.Abort()

    # report unrestrained energies, single point energy
    lmp.command("unfix REST")

    # omit energy contribution of the scanned entity
    #if omit_entity is True:
    #    null_entity(lmp, lmpdat, indices_and_values)

    try:
        lmp.command("run 0")
    except:
        print("***Error:  Simulation crashed (single step calculation)!")
        MPI.COMM_WORLD.Abort()

    lmp.close()


def measure_geometry(lmpdat, dcd, atm_ids):
    """
    Bla.

    Measure the actual bond, angle, dihedral achieved by restraining it.
    """
    sys_md = agum.Unification()
    sys_md.read_lmpdat(lmpdat)
    sys_md.import_dcd(dcd)
    sys_md.read_frames(frame=-2, to_frame=-1, frame_by="index")  # read only the last frame
    # decrement atom indices since for sys_md they start with 0, whereas
    # for lammps and gaussian they start with 1
    atm_idxs = [i - 1 for i in atm_ids]

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


def get_entities(gau_log, geom_entities):
    """
    Bla.

    Extract all entities from gaussian output file.
    """
    # pull all other entities (bonds, angles, dihedrals) from gaussian output
    sys_gau = agum.Unification()
    sys_gau.read_gau_log(gau_log)
    entities = []

    for geom_entity in geom_entities:
        cur_entities = sys_gau.gaussian_other_info[geom_entity]
        entities.append(cur_entities)

    return (entities, len(entities[0]))


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

    return "Angle [Degrees]"


def write_energies_file_header(filename):
    """
    Bla.

    Write header for energies file.
    """
    with open(filename, "w") as opened_filename:
        opened_filename.write("{:>20s} {:>20s}\n".format("Angle [Degrees]", "Energy [eV]"))


def guess_k(atom_ids):
    """
    Get k according to geometry entity.

    Parameters
    ----------
    atom_ids : list
        the ids of the atoms forming an entity, e.g. [25, 26, 28, 29] for a dihedral

    Returns
    -------
    k : float
        force constant k according to geometry type (bond, angle or dihedral)

    # shift phase by 180 degrees as defined in lamps manual
    # See: https://lammps.sandia.gov/doc/fix_restrain.html, "dihedral"

    """
    k = None

    if len(atom_ids) == 2:
        k = 1200
    elif len(atom_ids) == 3:
        k = 600
    elif len(atom_ids) == 4:
        k = 80
    else:
        raise Warning("Geometry entry seems odd")

    return k


def extract_atm_ids(geometries):
    """
    Extract all ids of the geometry strings in geometries.

    Parameters
    ----------
    geometries : list of str
        list of all geometries to scan simultaneously

    Returns
    -------
    ids : list of list of int
        all ids of all geometries to restrain, e.g. [[25, 26, 28, 29] [25, 26, 28, 30]]

    """
    atm_ids = []
    for geometry in geometries:
        geometry_ids = [int(i) for i in re.findall(r"\d+", geometry)]
        atm_ids.append(geometry_ids)

    return atm_ids


def md_from_ab_initio(gau_log, lmpdat, scanned_geom=None, add_geoms=None, temps=None, force_constants=None, energy_file_out="defaultname", output_idx=0):
    """
    Calculate md energy.

    Extract all necessary info from a gaussian log file, perform a simulated
    annealing simulation to get a relaxed structure and perform a minimization
    afterwards. This routine follows the path of the lowest energy, keeping
    the bond, angle or dihedral of interest the closest to that outcome of the
    ab initio calculation. Increasing the temperature or in-/decreasing the
    force constant k may be necessary.

    See: https://lammps.sandia.gov/doc/fix_restrain.html

    gau_log     str;
        gaussian log file
    lmpdat      str;
        lammps data file
    output      str;
        output-prefix for the data being written
    temp        tuple;
        starting and finishing temperature for the annealing procedure
    k           tuple;
        starting and stopping force constants for the geometric entity being examined
    scanned_geom : str
        the entity to scan as it is named in the gaussian output file, e.g.
        'A(7,8,25)' or 'D(5,7,11,13)'
    add_geoms : list of str
        geometries to add to the scanned one in order to restrain both simultaneously

    """
    # write output file header
    if not os.path.isfile(energy_file_out):
        write_energies_file_header(energy_file_out)

    # grab scanned entity from ab initio output and involved atom ids
    if scanned_geom is None:
        scanned_geom, _ = get_scanned_geometry(gau_log)

    if temps is None:
        temps = [600.0, 0.0]

    # include additional desired geometries that were not ab initio scanned,
    # but which the user wants to be restrained as well
    geometries = []
    geometries.append(scanned_geom)

    if add_geoms is not None:
        geometries.extend(add_geoms)

    # extract the plain atom ids from the entity strings
    total_geom_ids = extract_atm_ids(geometries)

    # get scanned geom atom ids for file naming
    scanned_geom_atm_ids = "_".join([str(i) for i in total_geom_ids[0]])

    # make an educated guess for each force constant of each geometry
    force_constants_start = 0.0

    if force_constants is None:
        force_constants = []

        for geom_ids in total_geom_ids:
            force_constants_end = guess_k(geom_ids)
            force_constants.append([force_constants_start, force_constants_end])

    # get ab initio values for each entity and number of geometries to restrain,
    # i.e. number of tasks to do
    total_scan_entities, tasks = get_entities(gau_log, geometries)

    # shift entity value by 180 degrees for dihedrals
    for geom_ids, idx in zip(total_geom_ids, xrange(len(total_scan_entities))):
        if len(geom_ids) == 4:
            total_scan_entities[idx] = [round(i + 180, 6) for i in total_scan_entities[idx]]

    for task_id in xrange(tasks):
        # each rank does its task and skips tasks assigned by other ranks (this
        # makes it possible to do restrain md runs in parallel)
        if (task_id % size != rank) or (rank > tasks):
            continue

        # define a dictionary with atom ids as key and the current
        # geometry value as value
        cur_geom_atm_ids_geom_value = {}

        for idx, geom_ids in enumerate(total_geom_ids):
            geom_entity = total_scan_entities[idx][task_id]
            geom_ids_str = " ".join([str(i) for i in geom_ids])
            cur_geom_atm_ids_geom_value[geom_ids_str] = geom_entity

        # define individual file name endings for each task and input gaussian file
        output_appendix = "{}_{}_{}".format(scanned_geom_atm_ids, task_id, output_idx)

        # skip calculations that have already been carried out
        if os.path.isfile(output_appendix + ".lmplog"):
            print(output_appendix + ".lmplog" + " already exists!")
            continue

        # annealing with quenching and minimization using lammps
        scan(lmpdat, output_appendix, cur_geom_atm_ids_geom_value, force_constants, temps=temps)

        # since minimized md-structure != minimized ab initio structure,
        # the real value of the geometry of interest must be derived from the
        # last frame of the md simulation
        current_geom_val = measure_geometry(lmpdat, output_appendix + ".dcd", total_geom_ids[0])

        # pull energy from md run
        last_pe_value = get_last_pe_value(output_appendix + ".lmplog")

        # pre defined string for each row in results file
        resume_file_row = "{:> 20.4f} {:> 20.8f}\n"

        with open(energy_file_out, "a") as opened_resume_file:
            opened_resume_file.write(resume_file_row.format(current_geom_val,
                                                            last_pe_value))


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
                            nargs="+",
                            help="Gaussian log-file")

        parser.add_argument("-geom_entity",
                            default=None,
                            help="Atom indices, where each string is one set of one geometry")

        parser.add_argument("-add_geom_entities",
                            default=None,
                            nargs="*",
                            help="Geometry to scan by gaussian syntax, e.g. A(2,1,3) A(2,1,4)")

        parser.add_argument("-k",
                            type=float,
                            default=None,
                            nargs="*",
                            help="Force Konstant K in eV. If add_geom_entities was chosen, a force constant for each entity has to be supplied as well. ")

        parser.add_argument("-out",
                            default="DEFAULTNAME",
                            help="Name of energy files.")

        args = parser.parse_args()

        # split k to sublists with two elements in each sublist (needed to create the right strings)
        if args.k is not None:
            k = []

            for i in args.k[::2]:
                for j in args.k[1::2]:
                    l = [i, j]
                k.append(l)

            args.k = k
    else:
        args = None

    args = comm.bcast(args, 0)

    #==============================================================================#
    # do restrained md for geometry scanned in gaussian
    #==============================================================================#
    output_file = "{}_md.txt".format(args.out)
    #pdb.set_trace()

    if not os.path.isfile(output_file):
        for gau_file_idx, cur_gau_log in enumerate(args.gau_logs):
            md_from_ab_initio(cur_gau_log, args.lmpdat, scanned_geom=args.geom_entity, add_geoms=args.add_geom_entities, energy_file_out=output_file, output_idx=gau_file_idx, force_constants=args.k)

    # wait for all ranks to finish
    time.sleep(2)
    print(bcolors.blue + "{} is done".format(rank) + bcolors.endc)

    # norm energies if file does not exist yet
    normed_md_file = "{}_md_normed.txt".format(args.out)

    if rank == 0 and not os.path.isfile(normed_md_file):
        norm_energy(output_file, normed_md_file)
