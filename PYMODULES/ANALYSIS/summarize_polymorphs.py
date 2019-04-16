#!/usr/bin/env python

"""
TODO EACH ITERATION OR AB INITIO OR EXPERIMENT SHOULD GET ITS OWN CELL IN
TODO JUPYTER NOTEBOOK AND EXECUTED THERE
TODO THIS PREVENTS CONFUSION AND MAKES EVERYTHING EASIER TO UNDERSTAND
"""

from __future__ import print_function, division
import numpy as np
import pandas as pd
import pdb
import md_box as mdb
from ag_pw import read_pw_out
from top_crds2lmpdat import get_cell_from_cif, get_other_stuff_from_cif
import ag_lmplog
import ag_unify_md as agum


class Polymorph(object):
    """
    """
    def __init__(self, atoms_p_molecule):
        """
        """
        self.box = None
        self.energy = None
        self.density = None
        self.volume = None
        self.temp = None
        self.natoms = None
        self.energy_kcal = None
        self.ucell_a = None
        self.ucell_b = None
        self.ucell_c = None
        self.energy_p_molecule = None
        self.energy_gain = None
        self.volume_p_molecule = None
        self.nmols = None
        self.density_p_molecule = None
        self.atoms_p_molecule = atoms_p_molecule

    def read_file(self, filename, filetype):
        """Get values of interest from given file of given type."""
        if filetype == "pwout":
            molsys = read_pw_out(filename)
            molsys.ts_boxes[-1].box_cart2lat()
            # get box and convert angles to degrees
            self.box = molsys.ts_boxes[-1]
            self.box.ltc_alpha = np.degrees(molsys.ts_boxes[-1].ltc_alpha)
            self.box.ltc_beta = np.degrees(molsys.ts_boxes[-1].ltc_beta)
            self.box.ltc_gamma = np.degrees(molsys.ts_boxes[-1].ltc_gamma)

            self.energy = molsys.pw_other_info["ENERGIES"][-1]
            self.density = molsys.pw_other_info["DENSITIES"][-1]
            self.volume = molsys.pw_other_info["VOLUMES"][-1]
            # ab initio is always at 0 K
            self.temp = 0.0
            self.natoms = len(molsys.ts_coords[-1])
        elif filetype == "lmplog":
            molsys = ag_lmplog.LmpLog()
            molsys.read_lmplog(filename)

            # box stuff
            self.box = mdb.Box()
            self.box.ltc_a = molsys.data[-1]["Cella"][-1]
            self.box.ltc_b = molsys.data[-1]["Cellb"][-1]
            self.box.ltc_c = molsys.data[-1]["Cellc"][-1]
            self.box.ltc_alpha = molsys.data[-1]["CellAlpha"][-1]
            self.box.ltc_beta = molsys.data[-1]["CellBeta"][-1]
            self.box.ltc_gamma = molsys.data[-1]["CellGamma"][-1]

            # other
            self.energy = molsys.data[-1]["PotEng"][-1]
            self.density = molsys.data[-1]["Density"][-1]
            self.volume = molsys.data[-1]["Volume"][-1]
            self.temp = molsys.data[-1]["Temp"][-1]

            with open(filename) as fin:
                line = fin.readline()
                while line != '':
                    if line.startswith("Loop time of"):
                        self.natoms = int(line.split()[-2])

                    line = fin.readline()

        elif filetype == "cif":
            # box stuff
            self.box = get_cell_from_cif(filename)
            cif_others = get_other_stuff_from_cif(filename)
            self.box.ltc_alpha = np.degrees(self.box.ltc_alpha)
            self.box.ltc_beta = np.degrees(self.box.ltc_beta)
            self.box.ltc_gamma = np.degrees(self.box.ltc_gamma)

            # other
            self.density = cif_others["density"]
            self.volume = cif_others["cell_volume"]
            self.temp = cif_others["cell_measurement_temperature"]
        else:
            raise Warning("Only pwout, lmplog and cif are valid keywords")

        if filetype == "lmplog" or filetype == "pwout":
            self.nmols = int(self.natoms / self.atoms_p_molecule)
        else:
            self.nmols = cif_others["cell_formula_units"]

    def convert_energy(self, energy_unit="kcal_mol"):
        """
        Convert the energy from a given file to the energy unit given.
        """
        if energy_unit == "kcal_mol":
            self.energy *= 23

    def get_density_p_molecule(self):
        """
        Density per molecule.
        """
        self.density_p_molecule = self.density / self.nmols

    def get_energy_p_molecule(self):
        """
        Calculate the energy per molecule.
        """
        self.energy_p_molecule = self.energy / self.nmols

    def get_energy_generation(self, energy_pmol_in_vacuo):
        """
        Get the energy that the crystallization generates.
        """
        self.energy_gain = self.energy_p_molecule - energy_pmol_in_vacuo

    def get_volume_p_molecule(self):
        """
        Calculate the volume per molecule.
        """
        self.volume_p_molecule = self.volume / self.nmols

    def supercell_to_unitcell(self, unitcell_box, supercell_box=None):
        """
        Reduce the box vector side lengths to that of a unit cell.
        To get the factor one could also use the original unchanged
        supercell-box if the box stretched somewhat too strong in comparison
        to its original box.
        """
        if supercell_box is None:
            cbox = self.box
        else:
            cbox = supercell_box

        factor_a = int(round(cbox.ltc_a / unitcell_box.ltc_a))
        factor_b = int(round(cbox.ltc_b / unitcell_box.ltc_b))
        factor_c = int(round(cbox.ltc_c / unitcell_box.ltc_c))
        self.ucell_a = self.box.ltc_a / factor_a
        self.ucell_b = self.box.ltc_b / factor_b
        self.ucell_c = self.box.ltc_c / factor_c


if __name__ == '__main__':
    ##################
    #  experiment    #
    ##################
    main_path = "/home/gadelmeier/SSHFS/hades/Research.new/carbamazepine/1.experiment/{}"
    exp_cbzi   = main_path.format("polymorph_I/1.CBZI.cif")
    exp_cbzii  = main_path.format("polymorph_II/1.CBZII.cif")
    exp_cbziii = main_path.format("polymorph_III/1.CBZIII.cif")
    exp_cbziv  = main_path.format("polymorph_IV/1.CBZIV.cif")
    exp_cbzv   = main_path.format("polymorph_V/1.CBZV.cif")
    #cif_ucell_files  = (exp_cbzi, exp_cbzii, exp_cbziii, exp_cbziv, exp_cbzv)
    cif_ucell_files  = {"I": (exp_cbzi, 6.1), "II": (exp_cbzii, 5.72), "III": (exp_cbziii, 6.41), "IV": (exp_cbziv, 5.95), "V": (exp_cbzv, 0.0)}

    exp_polymorphs = {}
    # energies in kcal*mol-1 from GRZESIAK, ADAM L., et al.,
    # JOURNAL OF PHARMACEUTICAL SCIENCES, vol. 92, no. 11, 2003, p. 12.
    #cif_energies = (6.1, 5.72, 6.41, 5.95, 0.0)

    for polymorph, (ucell_file, energy) in cif_ucell_files.iteritems():
        cur_polymorph = Polymorph(30)
        cur_polymorph.read_file(ucell_file, "cif")
        cur_polymorph.get_volume_p_molecule()
        cur_polymorph.energy_gain = energy
        exp_polymorphs[polymorph] = cur_polymorph

    del cur_polymorph

    ####################################
    #  ab initio (quantum espresso)    #
    ####################################
    pw_cbz_single = "/home/gadelmeier/SSHFS/hades/Research.new/carbamazepine/2.ab_initio/1.geom_opt/1.single/2.quantum_espresso/CBZIII_residue0_box_40_40_40_gamma/CBZIII_residue0_box_40_40_40.pwscf_out"
    PW_SINGLE_MOLECULE = agum.Unification()
    PW_SINGLE_MOLECULE.read_pwout(pw_cbz_single)

    main_path = "/home/gadelmeier/SSHFS/hades/Research.new/carbamazepine/2.ab_initio/1.geom_opt/3.unit_cell/3.vc_relax/{}"
    pw_cbzi   = main_path.format("CBZI/CBZI_vc_relax.pwscf_out")
    pw_cbzii  = main_path.format("CBZII/CBZII_vc_relax.pwscf_out")
    pw_cbziii = main_path.format("CBZIII/CBZIII_vc_relax.pwscf_out")
    pw_cbziv  = main_path.format("CBZIV/CBZIV_vc_relax.pwscf_out")
    pw_cbzv   = main_path.format("CBZV/CBZV_vc_relax.pwscf_out")
    #pw_ucell_files  = (pw_cbzi, pw_cbzii, pw_cbziii, pw_cbziv, pw_cbzv)
    pw_ucell_files  = {"I": pw_cbzi, "II": pw_cbzii, "III": pw_cbziii, "IV": pw_cbziv, "V": pw_cbzv}

    abinitio_polymorphs = {}

    for polymorph, cfile in pw_ucell_files.iteritems():
        cur_polymorph = Polymorph(30)
        cur_polymorph.read_file(cfile, "pwout")
        cur_polymorph.get_volume_p_molecule()
        # convert to kcal*mol-1
        cur_polymorph.convert_energy()
        cur_polymorph.get_energy_p_molecule()
        energy_p_molecule_vacuo = PW_SINGLE_MOLECULE.pw_other_info["ENERGIES"][-1] * 23
        cur_polymorph.get_energy_generation(energy_p_molecule_vacuo)
        abinitio_polymorphs[polymorph] = cur_polymorph
        del cur_polymorph

    ##############################################################
    # molecular dynamics (with different force field iterations) #
    ##############################################################
    path_single = "/home/gadelmeier/SSHFS/hades/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/1.single_molecule/CBZ_gaff-{0}_annealing/CBZ_gaff-{0}_annealing.lmplog"

    lammps_iterations = {}
    main_path = "/home/gadelmeier/SSHFS/hades/Research.new/carbamazepine/3.1.force_field_gaff/2.geom_opt/3.unit_cell/3.cell_relax/{}"
    path_scell_dreiding_off = main_path.format("CBZ_gaff-{0}/CBZ{1}/min_dreiding_off/CBZ{1}_gaff-{0}-min.lmplog")
    path_scell_dreiding_on  = main_path.format("CBZ_gaff-{0}/CBZ{1}/min_dreiding_on/CBZ{1}_gaff-{0}-min.lmplog")
    lmpdat_supercell        = main_path.format("CBZ_gaff-{0}/CBZ{1}/CBZ{1}_gaff-{0}_supercell.lmpdat")

    iterations_on = (0, 107)
    #iterations_off = (0, 107, 115)
    dreiding = "off"

    for iteration in iterations_on:
        # get box from original (multiplied and untouched) supercell in order to get the
        # factors of each unit cell it was built from
        # using the box from the last snapshot could lead to inconsistencies due to rounding
        #original_supercell = agum.Unification()
        #original_supercell.read_lmpdat(lmpdat_supercell.format(iteration, polymorph))
        #original_supercell.ts_boxes[0].box_lmp2lat()

        lammps_polymorphs = {}

        FF_SINGLE_MOLECULE = ag_lmplog.LmpLog()
        FF_SINGLE_MOLECULE.read_lmplog(path_single.format(iteration))
        # convert to kcal*mol-1
        FF_SINGLE_MOLECULE.data[-1]["PotEng"][-1] *= 23

        for index, polymorph in enumerate(("I", "II", "III", "IV", "V")):

            if iteration == 0 or dreiding == "off":
                # no dreiding force field for stock gaff force field
                cur_file = path_scell_dreiding_off.format(iteration, polymorph)
            else:
                cur_file = path_scell_dreiding_on.format(iteration, polymorph)

            #print(cur_file)
            cur_polymorph = Polymorph(30)
            cur_polymorph.read_file(cur_file, "lmplog")
            # convert to kcal*mol-1
            cur_polymorph.convert_energy()
            cur_polymorph.get_energy_p_molecule()
            cur_polymorph.get_volume_p_molecule()
            cur_polymorph.get_energy_generation(FF_SINGLE_MOLECULE.data[-1]["PotEng"][-1])
            cur_polymorph.supercell_to_unitcell(exp_polymorphs[polymorph].box)
            lammps_polymorphs[polymorph] = cur_polymorph

        lammps_iterations[iteration] = lammps_polymorphs

    # create pandas table
    rows = ["", "a", "b", "c", "alpha", "beta", "gamma", "vol p. molecule / A**3",
            "density / g*cm**-3", "E lattice p. molecule / kcal*mol**-1", "Energy gain / kcal*mol**-1"]
    df_combined = pd.DataFrame()
    #dataframes = {}

    #frames = []
    frames = {}
    frames["exp"] = []
    frames["abi"] = []
    frames["ff"]  = {}
    columnwidth = "{}"

    for polymorph in ("I", "II", "III", "IV", "V"):
        #frames = []

        exp_dataframe = pd.DataFrame(
            [
                "experiment",
                exp_polymorphs[polymorph].box.ltc_a,
                exp_polymorphs[polymorph].box.ltc_b,
                exp_polymorphs[polymorph].box.ltc_c,
                exp_polymorphs[polymorph].box.ltc_alpha,
                exp_polymorphs[polymorph].box.ltc_beta,
                exp_polymorphs[polymorph].box.ltc_gamma,
                exp_polymorphs[polymorph].volume_p_molecule,
                exp_polymorphs[polymorph].density,
                exp_polymorphs[polymorph].energy,
                exp_polymorphs[polymorph].energy_gain,
            ], index=rows, columns=[polymorph])

        frames["exp"].append(exp_dataframe)

        abinitio_dataframe = pd.DataFrame(
            [
                "ab initio PBE (USPP/Meyer)",
                abinitio_polymorphs[polymorph].box.ltc_a,
                abinitio_polymorphs[polymorph].box.ltc_b,
                abinitio_polymorphs[polymorph].box.ltc_c,
                abinitio_polymorphs[polymorph].box.ltc_alpha,
                abinitio_polymorphs[polymorph].box.ltc_beta,
                abinitio_polymorphs[polymorph].box.ltc_gamma,
                abinitio_polymorphs[polymorph].volume_p_molecule,
                abinitio_polymorphs[polymorph].density,
                abinitio_polymorphs[polymorph].energy_p_molecule,
                abinitio_polymorphs[polymorph].energy_gain,
            ], index=rows, columns=[polymorph])

        frames["abi"].append(abinitio_dataframe)
        #pdb.set_trace()

        for iteration in lammps_iterations:

            # create new list if key does not exist already
            try:
                frames["ff"][iteration]
            except KeyError:
                frames["ff"][iteration] = []

            forcefield_dataframe = pd.DataFrame(
                [
                    "gaff-{} dreiding {}".format(iteration, dreiding),
                    lammps_iterations[iteration][polymorph].ucell_a,
                    lammps_iterations[iteration][polymorph].ucell_b,
                    lammps_iterations[iteration][polymorph].ucell_c,
                    lammps_iterations[iteration][polymorph].box.ltc_alpha,
                    lammps_iterations[iteration][polymorph].box.ltc_beta,
                    lammps_iterations[iteration][polymorph].box.ltc_gamma,
                    lammps_iterations[iteration][polymorph].volume_p_molecule,
                    lammps_iterations[iteration][polymorph].density,
                    lammps_iterations[iteration][polymorph].energy_p_molecule,
                    lammps_iterations[iteration][polymorph].energy_gain,
                ], index=rows, columns=[polymorph])

            frames["ff"][iteration].append(forcefield_dataframe)

        #break

    table_exp = pd.concat(frames["exp"], axis=1)
    table_abi = pd.concat(frames["abi"], axis=1)
    pd.set_option('display.max_colwidth', 20)
    display(table_exp)
    display(table_abi)

    for key, iteration in frames["ff"].iteritems():
        ctable = pd.concat(iteration, axis=1)
        #print(ctable)
        display(ctable)
