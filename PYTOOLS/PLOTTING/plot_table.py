#!/usr/bin/env python

import os
import re
import argparse
import numpy as np
import pandas as pd
import pdb
from ag_pw import read_pw_out
from top_crds2lmpdat import get_cell_from_cif, get_other_stuff_from_cif
import ag_lmplog
from collections import OrderedDict


def get_values(filename, filetype):
    """Get values of interest from given file of given type."""
    box_a = 0.0
    box_b = 0.0
    box_c = 0.0
    alpha = 0.0
    beta = 0.0
    gamma = 0.0
    energy = 0.0
    density = 0.0
    volume = 0.0
    temp = 0.0

    if filetype == "pwout":
        molsys = read_pw_out(filename)
        molsys.ts_boxes[-1].box_cart2lat()
        box = molsys.ts_boxes[-1]
        box_a = box.ltc_a
        box_b = box.ltc_b
        box_c = box.ltc_c
        alpha = np.degrees(box.ltc_alpha)
        beta = np.degrees(box.ltc_beta)
        gamma = np.degrees(box.ltc_gamma)
        energy = molsys.pw_other_info["ENERGIES"][-1]
        density = molsys.pw_other_info["DENSITIES"][-1]
        volume = molsys.pw_other_info["VOLUMES"][-1]
        temp = 0.0  # ab initio is always at 0 K
    elif filetype == "lmplog":
        molsys = ag_lmplog.LmpLog()
        molsys.read_lmplog(filename)
        box_a = molsys.data[-1]["Cella"][-1]
        box_b = molsys.data[-1]["Cellb"][-1]
        box_c = molsys.data[-1]["Cellc"][-1]
        alpha = molsys.data[-1]["CellAlpha"][-1]
        beta = molsys.data[-1]["CellBeta"][-1]
        gamma = molsys.data[-1]["CellGamma"][-1]
        energy = molsys.data[-1]["PotEng"][-1]
        density = molsys.data[-1]["Density"][-1]
        volume = molsys.data[-1]["Volume"][-1]
        temp = molsys.data[-1]["Temp"][-1]
    elif filetype == "cif":
        box = get_cell_from_cif(filename)
        cif_others = get_other_stuff_from_cif(filename)
        box_a = box.ltc_a
        box_b = box.ltc_b
        box_c = box.ltc_c
        alpha = np.degrees(box.ltc_alpha)
        beta = np.degrees(box.ltc_beta)
        gamma = np.degrees(box.ltc_gamma)
        density = cif_others["density"]
        volume = cif_others["cell_volume"]
        temp = cif_others["cell_measurement_temperature"]
    else:
        raise Warning("Only pwout, lmplog and cif are valid keywords")

    keys = (
        "a",
        "b",
        "c",
        "alpha",
        "beta",
        "gamma",
        "temp",
        "volume",
        "density",
        "energy",
    )
    values = (box_a, box_b, box_c, alpha, beta, gamma, temp, volume, density, energy)
    return OrderedDict(list(zip(keys, values)))


def norm_values(keys_and_values, nmols, abc=(1, 1, 1)):
    """
    Transform original values where necessary.

    Parameters
    ----------
    keys_and_values : dict
        cell information
    nmols : int
        number of molecules to divide by
    abc : tuple of ints
        numbers to divide the super cell vectors a, b and c by

    """
    keys_and_values["a"] /= abc[0]
    keys_and_values["b"] /= abc[1]
    keys_and_values["c"] /= abc[2]
    keys_and_values["volume"] /= nmols
    keys_and_values["energy"] /= nmols
    return keys_and_values


def norm_by_polymorph(polymorph, keys_and_values, supercell=False):
    """Norm values according to the polymorph given."""
    # number of molecules for each unit cell
    nmols_p_unit_cell = {"I": 8, "II": 18, "III": 4, "IV": 8, "V": 8}
    # multiplication ratios of the unit cell as n-times a, b, c
    unit_cells_per_supercell = {
        "I": [12, 3, 3],
        "II": [2, 2, 13],
        "III": [8, 5, 4],
        "IV": [2, 8, 4],
        "V": [6, 5, 2],
    }

    if supercell is True:
        nmols = (
            unit_cells_per_supercell[polymorph][0]
            * unit_cells_per_supercell[polymorph][1]
            * unit_cells_per_supercell[polymorph][2]
            * nmols_p_unit_cell[polymorph]
        )
        abc = unit_cells_per_supercell[polymorph]
    else:
        nmols = nmols_p_unit_cell[polymorph]
        abc = [1, 1, 1]

    normed_keys_and_values = norm_values(keys_and_values, nmols, abc)
    return normed_keys_and_values


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("polymorph")
    PARSER.add_argument("-cif", default=None)
    PARSER.add_argument("-lmplogs", nargs="*", default=None)
    PARSER.add_argument("-pw", default=None)
    PARSER.add_argument("-no_normation", action="store_true")
    # PARSER.add_argument("-iteration", default=None)
    ARGS = PARSER.parse_args()

    COLUMNS = []
    LIST_OF_SERIES = []

    # extract values
    if ARGS.cif is not None:
        COLUMNS.append("exp.")

        if ARGS.no_normation is True:
            RESULTS = get_values(ARGS.cif, "cif")
        else:
            RESULTS = norm_by_polymorph(ARGS.polymorph, get_values(ARGS.cif, "cif"))

        SERIES = pd.Series(RESULTS, index=list(RESULTS.keys()))
        LIST_OF_SERIES.append(SERIES)
        print(SERIES)

    if ARGS.pw is not None:
        COLUMNS.append("ab initio")

        if ARGS.no_normation is True:
            RESULTS = get_values(ARGS.pw, "pwout")
        else:
            RESULTS = norm_by_polymorph(ARGS.polymorph, get_values(ARGS.pw, "pwout"))

        SERIES = pd.Series(RESULTS, index=list(RESULTS.keys()))
        LIST_OF_SERIES.append(SERIES)
        print(SERIES)

    if ARGS.lmplogs is not None:
        for lmplog in ARGS.lmplogs:
            ITERATION = re.findall(r"\d+", os.path.basename(lmplog))[0]
            # COLUMNS.append("gaff-{}".format(ITERATION))
            COLUMNS.append("{}".format(os.path.basename(lmplog).rstrip(".lmplog")))

            if ARGS.no_normation is True:
                RESULTS = get_values(lmplog, "lmplog")
            else:
                RESULTS = norm_by_polymorph(
                    ARGS.polymorph, get_values(lmplog, "lmplog"), supercell=True
                )

                SERIES = pd.Series(RESULTS, index=list(RESULTS.keys()))

            print(SERIES)
            LIST_OF_SERIES.append(SERIES)

    #
    pd.options.display.float_format = "{:,.4f}".format
    df = pd.DataFrame(LIST_OF_SERIES, index=COLUMNS).transpose()
    print(df)
    df.to_csv("polymorph_{}.csv".format(ARGS.polymorph), sep=",")
