#!/home/angad/Programs/anaconda/bin/python

from __future__ import print_function, division
import sys
import os
home = os.getenv("HOME")
sys.path.append(home + "Python/myPYMODULES/DLPOPLY_MODULES/")
sys.path.append(home + "Python/myPYMODULES/OTHER_MODULES/")
sys.path.append(home + "Python/myPYMODULES/MATH_MODULES/")
import tqdm
import numpy as np
import DL_HISTORY_class8_new_approaches_8_2 as dlhc
import my_lin_alg_compendium as mlc
import collections

# Usage: script, HISTORY, FIELD, EVERY-NTH-RESIDUE (0=all), START-STEP (Optional)
script = sys.argv[0]
historyfile = sys.argv[1]
fieldfile = sys.argv[2]
residue_number = sys.argv[3]

# every residue: 0; allowed residues are: 1, 2, 3, 4
residue_number = int(residue_number)
residue_number_1 = int(residue_number)

# Time step (if given) from which on to start evaluation
try:
    step_start = sys.argv[4]
    step_start = int(step_start)
except IndexError:
    step_start = None  # process every timestep

write_stddev = False

# FIELD-INSTANCE ///////////////////////////////////////////////////////////////
finfo = dlhc.FieldFile(fieldfile)
finfo.read_field()

# HISTORY-TO-EVALUATE //////////////////////////////////////////////////////////
# get number of time-steps
hist = dlhc.HistoryFile(finfo, historyfile)
hist.read_nsteps(ts_stop=None)
process_keys = hist.nstep_keys[step_start:]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# prepare necessary data-containers for the several molecule-types and molecules
raw_coms = {}  # delta com-values
raw_chis = {}  # delta chi-values
raw_phis = {}  # delta phi-values
raw_psis = {}  # delta psi-values

# delta values for each molecule
delta_com_all_mols = {}
delta_chi_all_mols = {}
delta_phi_all_mols = {}
delta_psi_all_mols = {}
# stddevs for each molecule
stddev_com_all_mols = {}
stddev_chi_all_mols = {}
stddev_phi_all_mols = {}
stddev_psi_all_mols = {}

com_squared_dists = {}  # values of current com minus the mean of all coms (for current molecule)
chi_squared_dists = {}  # values of current chi minus the mean of all chis (for current molecule)
phi_squared_dists = {}
psi_squared_dists = {}

# one container for each molecule in the FIELD
for mol_type in finfo.molecule_types:
    # containers for raw data
    raw_coms[mol_type.molname] = []
    raw_chis[mol_type.molname] = []
    raw_phis[mol_type.molname] = []
    raw_psis[mol_type.molname] = []
    # containers for mean values of each molecule
    delta_com_all_mols[mol_type.molname] = []
    delta_chi_all_mols[mol_type.molname] = []
    delta_phi_all_mols[mol_type.molname] = []
    delta_psi_all_mols[mol_type.molname] = []
    # containers for cur_value - mean_value for each molecule
    com_squared_dists[mol_type.molname] = []
    chi_squared_dists[mol_type.molname] = []
    phi_squared_dists[mol_type.molname] = []
    psi_squared_dists[mol_type.molname] = []

    # containers for stddevs
    if write_stddev:
        stddev_com_all_mols[mol_type.molname] = {}
        stddev_chi_all_mols[mol_type.molname] = {}
        stddev_phi_all_mols[mol_type.molname] = {}
        stddev_psi_all_mols[mol_type.molname] = {}

    # one sub-container for each molecule
    for molec in xrange(mol_type.nummols):
        raw_coms[mol_type.molname].append([])
        raw_chis[mol_type.molname].append([])
        raw_phis[mol_type.molname].append([])
        raw_psis[mol_type.molname].append([])
        delta_com_all_mols[mol_type.molname].append([])
        delta_chi_all_mols[mol_type.molname].append([])
        delta_phi_all_mols[mol_type.molname].append([])
        delta_psi_all_mols[mol_type.molname].append([])
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# File-progress visualization stuff
pbar = tqdm.tqdm(process_keys)

# Generate Data for each molecule ++++++++++++++++++++++++++++++++++++++++++++++
# arbitrary vector to form angles chi, psi, phi
arb_ref = [1, 0, 0]
opened_hist = open(historyfile)  # HISTORY must be unwrapped before!

for cur_tstep, tqdmstep in zip(process_keys, pbar):
    # reset number for every timestep
    residue_number_1 = residue_number
    # read file...
    cur_hstep = dlhc.HistoryFile(finfo, opened_hist)
    # ...but only one step at a time
    cur_hstep.read_history(ts_start=cur_tstep, ts_stop=cur_tstep)

    for cur_hist_moltype in cur_hstep.timesteps[cur_tstep].config.molecule_types:

        # select molecule types and calculate coms and mcs
        if cur_hist_moltype.molecules:

            if "CBZ" in cur_hist_moltype.molname:
                cur_hist_moltype.get_all_coms()  # Center of masses
                cur_hist_moltype.get_all_mcs(25, 26, 13, 11)  # Int. Coord-sys
            elif "Something else" in cur_hist_moltype.molname:  # TBD
                cur_hist_moltype.get_all_coms()
            else:
                pass

        for cur_hist_molecule in cur_hist_moltype.molecules:

            # pick the corresponding molecules
            if residue_number_1 == 0:  # every residue
                pass
            elif (cur_hist_molecule.mol_id+1) % residue_number_1 == 0:
                # only every 4th residue, starting from residue_number
                residue_number_1 += 4
            else:
                # skip everything in between
                continue

            # convert type from "object" to "float64" so we can work with it
            if cur_hist_molecule.center_of_mass is not None:
                cur_hist_molecule.center_of_mass = cur_hist_molecule.center_of_mass.astype("float64")
                # stash unprocessed com data
                raw_coms[cur_hist_moltype.molname][cur_hist_molecule.mol_id].append(
                    cur_hist_molecule.center_of_mass)

            if cur_hist_molecule.coord_sys is not None:
                cur_hist_molecule.coord_sys = [i.astype("float64") for i in
                                               cur_hist_molecule.coord_sys]
                # define angles
                cur_chi = mlc.angle_between(cur_hist_molecule.coord_sys[0],
                                            arb_ref, deg=True)
                cur_phi = mlc.angle_between(cur_hist_molecule.coord_sys[1],
                                            arb_ref, deg=True)
                cur_psi = mlc.angle_between(cur_hist_molecule.coord_sys[2],
                                            arb_ref, deg=True)

                # debugging stuff
                #print("Current cur_phi: ", cur_phi)
                #print("Current chi:", cur_chi)
                #print("Current psi: ", cur_psi)

                # stash unprocessed angle data
                raw_chis[cur_hist_moltype.molname][cur_hist_molecule.mol_id].append(
                    cur_chi)
                raw_phis[cur_hist_moltype.molname][cur_hist_molecule.mol_id].append(
                    cur_phi)
                raw_psis[cur_hist_moltype.molname][cur_hist_molecule.mol_id].append(
                    cur_psi)

    # rewind file, since we have already passed the next line ("timestep line")
    # so we do not have to re-read everything from the start
    cur_hist_pos = opened_hist.tell()
    cur_hist_pos -= 100000  # seems to be enough to rewind
    opened_hist.seek(cur_hist_pos)
    # Progress visualization stuff
    pbar.set_description("Processing %s" % cur_tstep)
    pbar.update(1)

opened_hist.close()

# Make coms and angles of each molecule comparable to each other +++++++++++++++
# coms -------------------------------------------------------------------------
for mol_type in finfo.molecule_types:
    mol_id = 0  # reset mol_id for each molecule type

    # each molecule has timestep times center of masses (molecule_coms)
    for molecule_coms in raw_coms[mol_type.molname]:

        if not molecule_coms:  # empty list, skip
            continue

        # mean xyz-position of center of mass
        cur_mol_mean_com = np.mean(molecule_coms, axis=0, dtype='float64')

        # calculate distance for mean center of mass and current center of mass
        for cur_ts_com in molecule_coms:
            cur_molecule_squared = (
                np.linalg.norm(cur_ts_com - cur_mol_mean_com)
            )  # **2
            # stash delta distances
            com_squared_dists[mol_type.molname].append(cur_molecule_squared)
            delta_com_all_mols[mol_type.molname][mol_id].append(cur_molecule_squared)

        mol_id += 1

del mol_id

# angles -----------------------------------------------------------------------
for mol_type in finfo.molecule_types:
    mol_id = 0
    # each molecule has timestep-times angles
    for molecule_chis, molecule_phis, molecule_psis in zip(
            raw_chis[mol_type.molname],
            raw_phis[mol_type.molname],
            raw_psis[mol_type.molname]):

        if not molecule_chis:
            continue

        cur_mol_mean_chi = np.mean(molecule_chis)
        cur_mol_mean_phi = np.mean(molecule_phis)
        cur_mol_mean_psi = np.mean(molecule_psis)
        for cur_ts_chi, cur_ts_phi, cur_ts_psi in zip(
                molecule_chis,
                molecule_phis,
                molecule_psis):
            cur_mol_squared_chi = (cur_ts_chi - cur_mol_mean_chi)  # **2
            cur_mol_squared_phi = (cur_ts_phi - cur_mol_mean_phi)  # **2
            cur_mol_squared_psi = (cur_ts_psi - cur_mol_mean_psi)  # **2
            chi_squared_dists[mol_type.molname].append(cur_mol_squared_chi)
            phi_squared_dists[mol_type.molname].append(cur_mol_squared_phi)
            psi_squared_dists[mol_type.molname].append(cur_mol_squared_psi)
            delta_chi_all_mols[mol_type.molname][mol_id].append(cur_mol_squared_chi)
            delta_phi_all_mols[mol_type.molname][mol_id].append(cur_mol_squared_phi)
            delta_psi_all_mols[mol_type.molname][mol_id].append(cur_mol_squared_psi)

        mol_id += 1

del mol_id

if write_stddev:
    # Calculate stddevs. for each molecule of each molecule type +++++++++++++++++++
    # iterate through each molecule type
    for cur_moltype in finfo.molecule_types:
        # iterate through each molecule
        for cur_mol_id in xrange(cur_moltype.nummols):

            # skip empty lists
            if not delta_com_all_mols[cur_moltype.molname][cur_mol_id]:
                continue
            else:
                pass

            # get stddevs.
            stddev_cur_coms = np.std(delta_com_all_mols[cur_moltype.molname][cur_mol_id])
            stddev_cur_chis = np.std(delta_chi_all_mols[cur_moltype.molname][cur_mol_id])
            stddev_cur_phis = np.std(delta_phi_all_mols[cur_moltype.molname][cur_mol_id])
            stddev_cur_psis = np.std(delta_psi_all_mols[cur_moltype.molname][cur_mol_id])
            # save all stddevs.
            stddev_com_all_mols[cur_moltype.molname][cur_mol_id] = stddev_cur_coms
            stddev_chi_all_mols[cur_moltype.molname][cur_mol_id] = stddev_cur_chis
            stddev_phi_all_mols[cur_moltype.molname][cur_mol_id] = stddev_cur_phis
            stddev_psi_all_mols[cur_moltype.molname][cur_mol_id] = stddev_cur_psis

    # Sort by stddevs.
    for cur_moltype in finfo.molecule_types:

        # Skip empty dicts
        if not stddev_com_all_mols[cur_moltype.molname]:
            continue

        # sort dictionaries by values
        stddevs_coms_sorted = collections.OrderedDict()
        stddevs_coms_sorted[cur_moltype.molname] = collections.OrderedDict(
            sorted(
                stddev_com_all_mols[cur_moltype.molname].items(),
                key=lambda t: t[1]
            )
        )

        stddevs_chis_sorted = collections.OrderedDict()
        stddevs_chis_sorted[cur_moltype.molname] = collections.OrderedDict(
            sorted(
                stddev_chi_all_mols[cur_moltype.molname].items(),
                key=lambda t: t[1]
            )
        )

        stddevs_phis_sorted = collections.OrderedDict()
        stddevs_phis_sorted[cur_moltype.molname] = collections.OrderedDict(
            sorted(
                stddev_phi_all_mols[cur_moltype.molname].items(),
                key=lambda t: t[1]
            )
        )

        stddevs_psis_sorted = collections.OrderedDict()
        stddevs_psis_sorted[cur_moltype.molname] = collections.OrderedDict(
            sorted(
                stddev_psi_all_mols[cur_moltype.molname].items(),
                key=lambda t: t[1]
            )
        )

    del (stddev_com_all_mols,
         stddev_chi_all_mols,
         stddev_phi_all_mols,
         stddev_psi_all_mols)

    # Write sorted std.devs. values to file
    with open("RESULTS_STDDEV.txt", "w") as stddev_file:
        for cur_moltype in stddevs_coms_sorted:

            stddev_file.write(
                "{0:s}\n".format(cur_moltype)
            )

            stddev_file.write(
                "{0:>6s}{1:>25s}{0:>8s}{2:>25s}{0:>8s}{3:>25s}{0:>8s}{4:>25s}\n".format(
                    "MOL-ID",
                    "COM-STD.DEV.",
                    "CHI-STD.DEV.",
                    "PHI-STD.DEV.",
                    "PSI-STD.DEV."
                )
            )

            for (com_mol_id,
                 chi_mol_id,
                 phi_mol_id,
                 psi_mol_id) in zip(stddevs_coms_sorted[cur_moltype],
                                    stddevs_chis_sorted[cur_moltype],
                                    stddevs_phis_sorted[cur_moltype],
                                    stddevs_psis_sorted[cur_moltype]):

                stddev_file.write(
                    "{4:>7d}{0: >25.15E}{5:>8d}{1: >25.15E}{6:>8d}{2: >25.15E}{7:>8d}{3: >25.15E}\n".format(
                        stddevs_coms_sorted[cur_moltype][com_mol_id],
                        stddevs_chis_sorted[cur_moltype][chi_mol_id],
                        stddevs_phis_sorted[cur_moltype][phi_mol_id],
                        stddevs_psis_sorted[cur_moltype][psi_mol_id],
                        com_mol_id,
                        chi_mol_id,
                        phi_mol_id,
                        psi_mol_id
                    )
                )

# RESULTS TO FILE //////////////////////////////////////////////////////////////
# Filename stuff +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if residue_number == 0:
    residue_filename = "all"
else:
    residue_filename = residue_number

# Write data +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Write data for each molecule type without distinguishing single molecules
with open("RESULTS_residue_%s.txt" % (residue_filename), "w") as result_file:
    for mol_type in finfo.molecule_types:
        result_file.write("%s\n" % mol_type.molname)
        result_file.write(
            "{0:>25s}{1:>25s}{2:>25s}{3:>25s}\n".format(
                "COM-SHIFT",
                "CHI-SHIFT",
                "PHI-SHIFT",
                "PSI-SHIFT"
            )
        )

        for cur_squared_com, cur_squared_chi, cur_squared_phi, cur_squared_psi in zip(
                com_squared_dists[mol_type.molname],
                chi_squared_dists[mol_type.molname],
                phi_squared_dists[mol_type.molname],
                psi_squared_dists[mol_type.molname]):
            result_file.write(
                "{0: >25.15E}{1: >25.15E}{2: >25.15E}{3: >25.15E}\n".format(
                    cur_squared_com,
                    cur_squared_chi,
                    cur_squared_phi,
                    cur_squared_psi,
                )
            )
