from __future__ import print_function, division
import pdb
import argparse
import numpy as np
from mpi4py import MPI
from lammps import lammps
from ag_lammps import read_lmpdat
import ag_unify_log as agul
from restrain_scan import norm_energy

#==============================================================================#
# Setup MPI
#==============================================================================#

comm = MPI.COMM_WORLD
size = comm.Get_size()  # number of processes in communicator
rank = comm.Get_rank()  # process' id(s) within a communicator
split = comm.Split(rank, key=0)


def get_shift_vector(lmpdat, dcd):
    """
    """
    dimer_sys = read_lmpdat(lmpdat, dcd)
    vt_N25_C26 = dimer_sys.ts_coords[-1][24]
    vt_C26_O27 = dimer_sys.ts_coords[-1][25] - dimer_sys.ts_coords[-1][27]
    vt_shift = np.cross(vt_N25_C26, vt_C26_O27)
    # normalize vt_shift
    return vt_shift / np.linalg.norm(vt_shift)


def scan_coordinates(lmprst, displace_atoms, frozen_atoms, vt_shift, output, save_step=100000, non_covalent=None):
    """
    """
    thermargs = ["step", "temp", "pe", "ebond", "eangle", "edihed", "eimp", "evdwl", "ecoul"]
    lmp = lammps(cmdargs=["-screen", "lmp_out.txt"])
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
    lmp.command("read_restart {}".format(lmprst))

    if non_covalent is not None:
        lmp.file(non_covalent)

    atoms = " ".join(map(str, displace_atoms))
    #pdb.set_trace()
    lmp.command("group displaced id {}".format(atoms))
    lmp.command("")

    if frozen_atoms is not None:
        lmp.command("group frozen id {}".format(" ".join(map(str, frozen_atoms))))
        lmp.command("fix freeze frozen setforce 0.0 0.0 0.0")

    lmp.command("thermo_style custom " + " ".join(thermargs))
    lmp.command("thermo_modify lost warn flush yes")
    lmp.command("thermo {}".format(save_step))
    lmp.command("fix MOMENT all momentum 100 linear 1 1 1 angular")
    lmp.command("dump DUMP all dcd {} {}.dcd".format(save_step, output))

    for i in xrange(-200, 400):
        if i < 0:
            vt_add = 0.05 * vt_shift
        else:
            vt_add = -0.05 * vt_shift
        lmp.command("displace_atoms displaced move {a[0]} {a[1]} {a[2]} units box".format(a=vt_add))

        try:
            lmp.command("minimize 1e-6 1e-9 2000000 100000")
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
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("lmpdat", help="Lammps' data-file with force field for cbz dimer system")
    PARSER.add_argument("lmprst", help="Lammps' restart file from minimization run of the cbz dimer system")
    PARSER.add_argument("dcd", help="dcd file from minimization run of the cbz dimer system")
    PARSER.add_argument("-non_covalent", help="file with pair coeffs for lammps and/or dreiding parameters")
    PARSER.add_argument("-o", default="DEFAULTNAME", help="Prefix of output files")
    ARGS = PARSER.parse_args()

    #minimize_folder = "/home/gadelmeier/SSHFS/hades/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/2.dimers/2H_anti_opt/0/"
    #lammps_data = minimize_folder + "CBZ_gaff-0_2H_dimer.lmpdat"
    #lammps_data_optimized_dcd = minimize_folder + "CBZ_gaff-0_2H_dimer.lmpdat_out.dcd"
    #lammps_restart = minimize_folder + "CBZ_gaff-0_2H_dimer.lmpdat_out.lmprst"
    SHIFT_VECTOR = get_shift_vector(ARGS.lmpdat, ARGS.dcd)
    #outname = "2H_flexible_scan"
    scan_coordinates(ARGS.lmprst, range(31, 61), [25, 26, 28, 55, 56, 58], SHIFT_VECTOR, ARGS.o, non_covalent=ARGS.non_covalent)
    write_summary(calculate_distances(ARGS.lmpdat, ARGS.o + ".dcd", [25], [55]), get_energies(ARGS.o + ".lmplog"))
    norm_energy("md_energies.txt", "md_energies_normed.txt")
