"""
Perform a relaxed molecule scan utilizing lammps.

This script is intended to push or pull a molecule or atom towards/away from
another molecule or atoms while freezing some atoms (relaxed scan). It is
intended to checkout and fit vdw-parameters to ab initio data (scans for ab initio
therefor have to be done the same way).
"""

from __future__ import print_function, division
import pdb
import os
import numpy as np
from mpi4py import MPI
from lammps import lammps
from ag_lammps import read_lmpdat
import ag_unify_md as agum
import ag_unify_log as agul
from restrain_scan import norm_energy

#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator
split = comm.Split(rank, key=0)


# check for dreiding-settings in settings file
def check_settings_file(settings_file):
    """Check if dreiding is in the settings file."""
    with open(settings_file) as fin:
        for line in fin:
            if line.startswith("pair_style") and "hbond/dreiding" in line:
                return True

    return False


def create_lmpdat(lmptop, xyz, output_name="foobar"):
    """
    Create the lammps data file for the dimer to scan from the dimer given.

    Parameters
    ----------
    lmptop : str
        lammps-data file of the monomer

    xyz : str
        xyz file of the system that shall be build from the monomer

    """
    sys_lmptop = agum.Unification()
    sys_lmptop.read_lmpdat(lmptop)
    sys_coords = agum.Unification()
    sys_coords.read_xyz(xyz)

    # replace original topology atom coordinates
    sys_lmptop.ts_coords = sys_coords.ts_coords

    # replicate topology if necessary
    num_topo_atms = len(sys_lmptop.atoms)
    num_topo_coords = len(sys_lmptop.ts_coords[-1])
    n = num_topo_coords / num_topo_atms

    if n.is_integer():
        # replicate topology
        sys_lmptop.add_topology_replicate(int(n - 1), refresh_bonds=True)
    else:
        raise Warning("Number of atoms in coordinates file is not a multiple of atoms in topology file!")

    sys_lmptop.change_indices(incr=1, mode="increase")
    sys_lmptop.write_lmpdat(output_name, cgcmm=True)


# shift vector
def get_shift_vector(lmpdat, atoms_cog1, atoms_cog2, dcd=None, frame_id=-1):
    """
    Calculate the vector for shifting both molecules.

    Parameters
    ----------
    lmpdat : str
        lammps data file

    dcd : str or None
        default: None
        dcd file from molecular dynamics run which supplies the latest
        coordinates

    atoms_cog1 : tuple or list
        atom indices that form the first center of geometry;
        if only one atom is given, it will be the center of geometry

    atoms_cog2 : tuple or list
        same as atoms_cog1 but for 2nd center of geometry

    Returns:
    --------
    vt_shift : np-array
        normalized shifting vector

    """
    dimer = read_lmpdat(lmpdat, dcd)

    if len(atoms_cog1) > 1:
        cog1 = dimer.get_cog(frame_id, *atoms_cog1)
        cog2 = dimer.get_cog(frame_id, *atoms_cog2)
    else:
        cog1 = dimer.ts_coords[frame_id][atoms_cog1[0]]
        cog2 = dimer.ts_coords[frame_id][atoms_cog2[0]]

    vt_shift = cog1 - cog2
    vt_shift /= np.linalg.norm(vt_shift)
    return vt_shift


def minimize(lmpdat, settings_file, lmpcfg=None, output="foobar"):
    thermargs = ["step", "temp", "pe", "ebond", "eangle", "edihed", "eimp", "evdwl", "ecoul"]

    lmp = lammps(cmdargs=["-screen", "lmp_out.txt"])
    lmp.command("log {}.lmplog".format(output))
    lmp.file(settings_file)

    # read file
    if os.path.isfile(lmpdat):
        lmp.command("read_data {}".format(lmpdat))
    else:
        raise Warning("Could not find data file.")

    if lmpcfg is not None:
        lmp.file(lmpcfg)

    # add dreiding parameters if a run with a dreiding containing hybrid
    # style is run
    if check_settings_file(settings_file) is True:
        lmp.command("variable E_hbond equal c_hb[2]")
        lmp.command("compute hb all pair hbond/dreiding/lj")
        lmp.command("variable n_hbond equal c_hb[1]")
        lmp.command("thermo_style custom " + " ".join(thermargs) + " v_n_hbond v_E_hbond")
    else:
        lmp.command("thermo_style custom " + " ".join(thermargs))

    lmp.command("thermo_modify lost warn flush yes")
    lmp.command("thermo {}".format(10))
    lmp.command("fix MOMENT all momentum 100 linear 1 1 1 angular")
    lmp.command("dump DUMP all dcd {} {}.dcd".format(10, output))

    # perform a minimization of the whole system
    #lmp.command("min_style {}".format(min_style))

    try:
        lmp.command("minimize 1e-6 1e-9 2000000 1000000")
    except:
        print("***Error:  Simulation crashed (minimization)!")
        MPI.COMM_WORLD.Abort()

    lmp.command("write_restart {}.lmprst".format(output))


def scan_coordinates(lmprst, displace_atoms, frozen_atoms, vt_shift, output,
                     settings_file, save_step=100000, lmpcfg=None):
    """
    Move one part of the molecule along a shifting vector with sp calculations.

    - Perform minimization run for the whole system using one of lammps' implemented
      minimization styles
    - Freeze a certain amount of atoms while moving one part of the system away
      from the rest while keeping the rest of the system flexible
    - After each movement step relax the non frozen part of the whole system

    Sources
    -------
    https://lammps.sandia.gov/doc/min_style.html

    Parameters
    ----------
    lmprst : str
        lammps restart file from a previous minimization run

    displace_atoms : tuple or list
        indices of atoms to be displaced

    frozen_atoms : tuple or list
        indices of atoms that should not be movable during the simulation

    output : str
        output name of the system

    save_step : int
        amount of steps after something will be written in the lammps-log file

    settings_file : str
        file with basic molecular dynamics settings; must be provided by the user

    lmpcfg : str
        file with force field parameters for the run, most often for vdw settings
    """
    thermargs = ["step", "temp", "pe", "ebond", "eangle", "edihed", "eimp", "evdwl", "ecoul"]

    lmp = lammps(cmdargs=["-screen", "lmp_out.txt"])
    lmp.command("log {}.lmplog".format(output))
    lmp.file(settings_file)

    # read file
    if os.path.isfile(lmprst):
        lmp.command("read_restart {}".format(lmprst))
    else:
        raise Warning("Could not find restart file.")

    if lmpcfg is not None:
        lmp.file(lmpcfg)

    # add dreiding parameters if a run with a dreiding containing hybrid
    # style is run
    if check_settings_file(settings_file) is True:
        lmp.command("variable E_hbond equal c_hb[2]")
        lmp.command("compute hb all pair hbond/dreiding/lj")
        lmp.command("variable n_hbond equal c_hb[1]")
        lmp.command("thermo_style custom " + " ".join(thermargs) + " v_n_hbond v_E_hbond")
    else:
        lmp.command("thermo_style custom " + " ".join(thermargs))

    lmp.command("thermo_modify lost warn flush yes")
    lmp.command("thermo {}".format(save_step))
    lmp.command("fix MOMENT all momentum 100 linear 1 1 1 angular")
    lmp.command("dump DUMP all dcd {} {}.dcd".format(save_step, output))

    # perform a minimization of the whole system
    try:
        lmp.command("minimize 1e-6 1e-9 2000000 1000000")
    except:
        print("***Error:  Simulation crashed (minimization)!")
        MPI.COMM_WORLD.Abort()

    # scanning of the atoms
    atoms = " ".join(map(str, displace_atoms))
    lmp.command("group displaced id {}".format(atoms))

    if frozen_atoms is not None:
        lmp.command("group frozen id {}".format(" ".join(map(str, frozen_atoms))))
        lmp.command("fix freeze frozen setforce 0.0 0.0 0.0")

    for i in xrange(1, 200):

        if i < 20:
            # push both parts together
            vt_add = 0.05 * vt_shift
        else:
            # push one part farther away
            vt_add = -0.05 * vt_shift

        lmp.command("displace_atoms displaced move {a[0]} {a[1]} {a[2]} units box".format(a=vt_add))

        try:
            lmp.command("minimize 1e-6 1e-9 2000000 1000000")
        except:
            print("***Error:  Simulation crashed (minimization)!")
            MPI.COMM_WORLD.Abort()


def calculate_distances(lmpdat, dcd, idxs_atm1, idxs_atm2):
    """
    Calculate the distance between the center of geometry of the two benzene rings.
    """
    distances = []
    dimer_sys = read_lmpdat(lmpdat, dcd, frame_idx_start=0)
    #pdb.set_trace()

    for frame_idx in xrange(0, len(dimer_sys.ts_coords)):

        if len(idxs_atm1) == 1:
            cog1 = dimer_sys.ts_coords[frame_idx][idxs_atm1[0]]
        else:
            cog1 = dimer_sys.get_cog(frame_idx, *idxs_atm1)

        if len(idxs_atm2) == 1:
            cog2 = dimer_sys.ts_coords[frame_idx][idxs_atm2[0]]
        else:
            cog2 = dimer_sys.get_cog(frame_idx, *idxs_atm2)

        distance = np.linalg.norm(cog1 - cog2)
        distances.append(distance)

    return distances


def get_energies(lmplog):
    """
    """
    pot_energies = []
    coul_energies = []
    vdw_energies = []
    intra_energies = []

    log_data = agul.LogUnification()
    log_data.read_lmplog(lmplog)

    for entry in log_data.data:
        pot_energies.append(entry["PotEng"][-1])
        coul_energies.append(entry["E_coul"][-1])
        vdw_energies.append(entry["E_vdwl"][-1])
        intra_energies.append(entry["E_bond"][-1] + entry["E_angle"][-1] +
                              entry["E_dihed"][-1] + entry["E_impro"][-1])

    return (pot_energies, coul_energies, vdw_energies, intra_energies)


def write_summary(distances, energies):
    """
    """
    row = "{:> 20.4f} {:> 20.8f} {:> 20.4f} {:> 20.4f} {:> 20.4f}\n"
    with open("md_energies.txt", "w") as opened_file:
        opened_file.write("{:>20s} {:>20s} {:>20s} {:>20s} {:>20s}\n".format(
            "Distance [Angstrom]",
            "PotEng [eV]",
            "E_coul [eV]",
            "E_vdwl [eV]",
            "E_intra [eV]"))
        for distance, pe, coule, vdwe, intrae in zip(distances, energies[0], energies[1], energies[2], energies[3]):
            opened_file.write(row.format(distance, pe, coule, vdwe, intrae))


if __name__ == "__main__":
    # force field and settings
    ITERATION = "111"
    DREIDING = "on"
    SETTINGS_FILE = "/home/gadelmeier/hdd/SSHFS/hades/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/md_settings/settings_dreiding_{}.lmpcfg".format(DREIDING)
    FF_FILE = "/home/gadelmeier/hdd/SSHFS/hades/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/md_settings/CBZ_gaff-{}_dreiding_{}.lmpcfg".format(ITERATION, DREIDING)
    LMPDAT = "/home/gadelmeier/hdd/SSHFS/hades/Research.new/carbamazepine/3.1.force_field_gaff/4.forcefields/CBZ_gaff-{}_novdw.lmpdat".format(ITERATION)

    # H0 dimer
    execute_0H_scan = True

    if execute_0H_scan is True:

        XYZ_H0 = "/home/gadelmeier/hdd/SSHFS/hades/Research.new/carbamazepine/2.ab_initio/1.geom_opt/2.dimers/0H_anti_opt/1.gaussian09/CBZ_0H_sp_wB97XD_cc-pVTZ.gau.out.xyz"
        # iteration specific force field and settings
        # scan of the H-Bonds of two Carbamazepine molecules
        EMIN_FOLDER = "/home/gadelmeier/hdd/SSHFS/hades/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/2.dimers/0H_anti_opt/gaff-{}_dreiding_{}/".format(ITERATION, DREIDING)
        EMIN_OUT = "CBZ_Dimer_anti_0H_gaff-{}_dreiding_{}_minimized".format(ITERATION, DREIDING)

        # MINIMIZATION
        LMPDAT_H0 = EMIN_FOLDER + EMIN_OUT + ".lmpdat"
        EMIN_DCD_H0 = EMIN_FOLDER + EMIN_OUT + ".dcd"
        EMIN_LMPRST_H0 = EMIN_FOLDER + EMIN_OUT + ".lmprst"
        print(EMIN_LMPRST_H0)

        try:
            os.mkdir(EMIN_FOLDER)
        except OSError:
            pass

        os.chdir(EMIN_FOLDER)
        create_lmpdat(LMPDAT, XYZ_H0, output_name=LMPDAT_H0)
        minimize(LMPDAT_H0, SETTINGS_FILE, FF_FILE, output=EMIN_OUT)

        WORKING_DIR = "/home/gadelmeier/hdd/SSHFS/hades/Research.new/carbamazepine/3.1.force_field_gaff/3.scans/2.dimer_scans/0H_anti_relaxed_scan/1.vertical_scan/CBZ_gaff-{}_dreiding_{}".format(ITERATION, DREIDING)

        # scan
        try:
            os.mkdir(WORKING_DIR)
        except OSError:
            pass

        os.chdir(WORKING_DIR)
        VT_SHIFT = get_shift_vector(LMPDAT_H0, [14, 15, 16, 18, 20, 22], [44, 45, 46, 48, 50, 52], dcd=EMIN_DCD_H0)
        scan_coordinates(
            EMIN_LMPRST_H0,
            range(31, 61), [14, 15, 16, 18, 20, 22, 44, 45, 46, 48, 50, 52],
            VT_SHIFT, "CBZ_Dimer_anti_0H_gaff-{}_dreiding_{}_scan".format(ITERATION, DREIDING),
            SETTINGS_FILE, lmpcfg=FF_FILE)

        write_summary(calculate_distances(LMPDAT_H0, "CBZ_Dimer_anti_0H_gaff-{}_dreiding_{}_scan.dcd".format(ITERATION, DREIDING), [14, 15, 16, 18, 20, 22], [44, 45, 46, 48, 50, 52]), get_energies("CBZ_Dimer_anti_0H_gaff-{}_dreiding_{}_scan.lmplog".format(ITERATION, DREIDING)))
        norm_energy("md_energies.txt", "md_energies_normed.txt")

    # H2 dimer
    execute_2H_scan = True

    if execute_2H_scan is True:
        XYZ_H2 = "/home/gadelmeier/hdd/SSHFS/hades/Research.new/carbamazepine/2.ab_initio/1.geom_opt/2.dimers/2H_anti_opt/1.gaussian09/CBZ_2H_sp_wB97XD_cc-pVTZ.gau.out.xyz"
        # scan of the H-Bonds of two Carbamazepine molecules
        EMIN_FOLDER = "/home/gadelmeier/hdd/SSHFS/hades/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/2.dimers/2H_anti_opt/gaff-{}_dreiding_{}/".format(ITERATION, DREIDING)
        EMIN_OUT = "CBZ_Dimer_anti_2H_gaff-{}_dreiding_{}_minimized".format(ITERATION, DREIDING)

        # MINIMIZATION
        LMPDAT_H2 = EMIN_FOLDER + EMIN_OUT + ".lmpdat"
        EMIN_DCD_H2 = EMIN_FOLDER + EMIN_OUT + ".dcd"
        EMIN_LMPRST_H2 = EMIN_FOLDER + EMIN_OUT + ".lmprst"

        try:
            os.mkdir(EMIN_FOLDER)
        except OSError:
            pass

        os.chdir(EMIN_FOLDER)
        create_lmpdat(LMPDAT, XYZ_H2, output_name=LMPDAT_H2)
        minimize(LMPDAT_H2, SETTINGS_FILE, FF_FILE, output=EMIN_OUT)

        WORKING_DIR = "/home/gadelmeier/hdd/SSHFS/hades/Research.new/carbamazepine/3.1.force_field_gaff/3.scans/2.dimer_scans/2H_anti_relaxed_scan/1.horizontal_scan/CBZ_gaff-{}_dreiding_{}".format(ITERATION, DREIDING)

        # scan
        try:
            os.mkdir(WORKING_DIR)
        except OSError:
            pass

        os.chdir(WORKING_DIR)
        VT_SHIFT = get_shift_vector(LMPDAT_H2, [25], [55], dcd=EMIN_DCD_H2)

        scan_coordinates(
            EMIN_LMPRST_H2,
            range(31, 61), [25, 26, 28, 55, 56, 58],
            VT_SHIFT, "CBZ_Dimer_anti_2H_gaff-{}_dreiding_{}_scan".format(ITERATION, DREIDING),
            SETTINGS_FILE, lmpcfg=FF_FILE)

        write_summary(
            calculate_distances(
                LMPDAT_H2, "CBZ_Dimer_anti_2H_gaff-{}_dreiding_{}_scan.dcd".format(ITERATION, DREIDING),
                [25], [55]),
            get_energies("CBZ_Dimer_anti_2H_gaff-{}_dreiding_{}_scan.lmplog".format(ITERATION, DREIDING)))

        norm_energy("md_energies.txt", "md_energies_normed.txt")
