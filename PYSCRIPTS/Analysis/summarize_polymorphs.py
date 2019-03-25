#!/usr/bin/env python
from __future__ import print_function, division
import os
import re
import argparse
import numpy as np
#import scipy.constants as sc
import pandas as pd
import pdb
import md_box as mdb
from ag_pw import read_pw_out
from top_crds2lmpdat import get_cell_from_cif, get_other_stuff_from_cif
import ag_lmplog
#from collections import OrderedDict


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
                while line != "":
                    if line.startswith("Loop time of"):
                        self.natoms = int(line.split()[-2])

        elif filetype == "cif":
            # box stuff
            self.box = get_cell_from_cif(filename)
            cif_others = get_other_stuff_from_cif(filename)
            self.box.ltc_alpha = np.degrees(self.box.ltc_alpha)
            self.box.ltc_beta = np.degrees(self.box.ltc_beta)
            self.box.ltc_gamma = np.degrees(self.box.ltc_gamma)

            # other
            self.box.density = cif_others["density"]
            self.box.volume = cif_others["cell_volume"]
            self.box.temp = cif_others["cell_measurement_temperature"]
        else:
            raise Warning("Only pwout, lmplog and cif are valid keywords")

        if filetype == "lmplog" or filetype == "pw":
            self.nmols = self.natoms / self.atoms_p_molecule
        else:
            self.nmols = cif_others["cell_formula_units"]

    def convert_energy(self, energy_unit="kcal_mol"):
        """
        Convert the energy from a given file to the energy unit given.
        """
        if energy_unit == "kcal_mol":
            self.energy_kcal = self.energy * 23

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
        self.energy_gain = self.energy - energy_pmol_in_vacuo

    def get_volume_p_molecule(self):
        """
        Calculate the volume per molecule.
        """
        self.volume_p_molecule = self.volume / self.nmols

    def supercell_to_unitcell(self, original_cell_box):
        """
        Reduce the box vector side lengths to that of a unit cell.
        """
        factor_a = int(round(self.box.ltc_a / original_cell_box.ltc_a))
        factor_b = int(round(self.box.ltc_b / original_cell_box.ltc_b))
        factor_c = int(round(self.box.ltc_c / original_cell_box.ltc_c))
        self.ucell_a = self.box.ltc_a / factor_a
        self.ucell_b = self.box.ltc_b / factor_b
        self.ucell_c = self.box.ltc_c / factor_c


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("cif", default=None)
    PARSER.add_argument("-lmplogs", nargs="2", default=None, help="First file is the main cell, second file is the single molecule in vacuo")
    PARSER.add_argument("-pws", nargs="2", default=None, help="First file is the main cell, second file is the single molecule in vacuo")
    ARGS = PARSER.parse_args()

    ATOMS_P_MOLECULE = 30
    CIF_POLYMORPH = Polymorph(ATOMS_P_MOLECULE)
    CIF_POLYMORPH.read_file(ARGS.cif, "cif")
    CIF_POLYMORPH.volume_p_molecule()

    CIF_RESULTS = {"a": CIF_POLYMORPH.box.ltc_a,
                   "b": CIF_POLYMORPH.box.ltc_b,
                   "c": CIF_POLYMORPH.box.ltc_c,
                   "alpha": CIF_POLYMORPH.box.ltc_alpha,
                   "beta": CIF_POLYMORPH.box.ltc_beta,
                   "gamma": CIF_POLYMORPH.box.ltc_gamma,
                   "volume p. molecule": CIF_POLYMORPH.volume_p_molecule}

    if ARGS.lmplogs[0] is not None and ARGS.lmplogs[1] is not None:
        LMPLOG_SINGLE_MOLECULE = Polymorph(1)
        LMPLOG_SINGLE_MOLECULE.convert_energy()

        LMPLOG_SUPERCELL = Polymorph(ATOMS_P_MOLECULE)
        LMPLOG_SUPERCELL.convert_energy()
        LMPLOG_SUPERCELL.get_energy_p_molecule()
        LMPLOG_SUPERCELL.get_volume_p_molecule()
        LMPLOG_SUPERCELL.get_energy_generation(LMPLOG_SINGLE_MOLECULE.energy)
        LMPLOG_SUPERCELL.supercell_to_unitcell(CIF_POLYMORPH.box)

        MD_RESULTS = {"a": LMPLOG_SUPERCELL.ucell_a,
                      "b": LMPLOG_SUPERCELL.ucell_b,
                      "c": LMPLOG_SUPERCELL.ucell_c,
                      "alpha": LMPLOG_SUPERCELL.box.ltc_alpha,
                      "beta": LMPLOG_SUPERCELL.box.ltc_beta,
                      "gamma": LMPLOG_SUPERCELL.box.ltc_gamma,
                      "volume p. molecule": LMPLOG_SUPERCELL.volume_p_molecule,
                      "energy gain p. molecule": LMPLOG_SUPERCELL.energy_gain}

    if ARGS.pws[0] is not None and ARGS.pws[1] is not None:
        PW_SINGLE_MOLECULE = Polymorph(1)
        PW_SINGLE_MOLECULE.convert_energy()

        PW_UNITCELL = Polymorph(ATOMS_P_MOLECULE)
        PW_UNITCELL.convert_energy()
        PW_UNITCELL.get_energy_p_molecule()
        PW_UNITCELL.get_volume_p_molecule()
        PW_UNITCELL.get_energy_generation(PW_SINGLE_MOLECULE.energy)
