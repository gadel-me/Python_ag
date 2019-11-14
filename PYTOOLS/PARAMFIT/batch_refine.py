#!/usr/bin/env python

import os
import sys
#import shutil as sl
import argparse
import subprocess32 as sp32
import time
import ag_amberparamfit as ampf
import ntpath

__version__ = "2017-08-11"

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="batch_refine_params.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Automatize the fitting process by processing each step in \
                 the user's given order. 'set_params.in/out files must be given \
                 for this to work for each parameter to be fit. \
                 Quantum-files should be weighted before the process.")

parser.add_argument("-prmtop",
                    required=True,
                    metavar="*.prmtop",
                    action="store",
                    help="First amber topology file to initiate the fitting process.",
                    )

parser.add_argument("-mol2",
                    required=True,
                    help="MOL2 for tleap. Atomic order must be the same as in prmtop."
                    )

parser.add_argument("-mdcrds",
                    required=True,
                    nargs="*",
                    metavar="*.mdcrd",
                    help="Amber coordinates files with conformations from \
                         relaxed molecule scan.",
                    )

parser.add_argument("-quantums",
                    required=True,
                    nargs="*",
                    metavar="*.dat",
                    help="Extracted single point energies from gaussian calculations. \
                          order must be the same as in amber coordinates file!",
                    )

parser.add_argument("-frcmod",
                    help="First frcmod file from antechamber output")

parser.add_argument("-forcefield",
                    default="gaff2",
                    help="Force Field to load into tleap to generate a new prmtop.")

args = parser.parse_args()

# define amber environment vars ------------------------------------------------
AMBERHOME   = os.getenv("AMBERHOME")
paramfit    = "{}/bin/paramfit".format(AMBERHOME)
antechamber = "{}/bin/antechamber".format(AMBERHOME)
tleap       = "{}/bin/tleap".format(AMBERHOME)
python      = sys.executable
workdir     = os.getcwd()

# absolute paths for all files
args.prmtop = os.path.abspath(args.prmtop)
args.mol2   = os.path.abspath(args.mol2)
args.frcmod = os.path.abspath(args.frcmod)
args.mdcrds = [os.path.abspath(i) for i in args.mdcrds]
args.quantums = [os.path.abspath(i) for i in args.quantums]

# get total number parameters to fit (needed for indices)
num_params = len(args.quantums)
cur_prmtop = args.prmtop  # current prmtop file (substituted after each run)

# paramfit ---------------------------------------------------------------------
frcmods = []
frcmods.append(args.frcmod)

for prm_idx in range(num_params):
    cur_quantum = args.quantums[prm_idx]
    cur_mdcrd   = args.mdcrds[prm_idx]

    # get number of structures that were fitted for current parameter
    # from quantum energy file
    prm_name = ntpath.basename(cur_quantum).rstrip("_energies.dat")

    with open(cur_quantum, "r") as f_in:
        nstructures = 0
        for line in f_in:
            # skip empty lines
            if line.strip() != "":
                nstructures += 1

    print("***Info: Current parameter to fit is {}".format(prm_name))
    # change to directory of parameter to fit and execute paramfit
    try:
        os.mkdir(prm_name)
    except OSError:
        pass

    os.chdir(prm_name)

    # run paramfit -------------------------------------------------------------
    prmfit = ampf.Paramfit(nstructures)

    # set parameters to fit
    set_params_in = "set_params.in"
    set_params_out = "set_params.out"

    if os.path.isfile(set_params_in) is False:
        prmfit.write_set_params_in(set_params_in)

    if os.path.isfile(set_params_out) is False:
        set_parms = sp32.call([paramfit,
                               "-i", set_params_in,
                               "-p", cur_prmtop])
        time.sleep(2)

    # get initial K value
    fit_k_in  = "fit_K.in"
    fit_k_out = "fit_K.out"

    if os.path.isfile(fit_k_in) is False:
        prmfit.write_fit_k(fit_k_in)

    if os.path.isfile(fit_k_out) is False:
        with open(fit_k_out, "w") as out:
            print("***Info: Executing paramfit (initial K)")
            fit_K = sp32.call([paramfit,
                               "-i", fit_k_in,
                               "-p", cur_prmtop,
                               "-c", cur_mdcrd,
                               "-q", cur_quantum],
                              stdout=out)
            if fit_K != 0:
                raise Warning("Something went wrong with paramfit during 'derivation of K'.")

    force_const_K = prmfit.read_fit_k(fit_k_out)

    # fit current parameter
    fit_params_in  = "fit_params.in"
    fit_params_out = "fit_params.out"
    fitted_energy  = "fit_output_energy.dat"

    if os.path.isfile(fit_params_in) is False:
        prmfit.write_fit_params(fit_params_in, set_params_out, force_const_K)

    if os.path.isfile(fit_params_out) is False:
        with open(fit_params_out, "w") as out:
            print("***Info: Executing paramfit (fitting parameters)")
            fit_param = sp32.call([paramfit,
                                   "-i", fit_params_in,
                                   "-p", cur_prmtop,
                                   "-c", cur_mdcrd,
                                   "-q", cur_quantum],
                                  stdout=out)

            if fit_param != 0:
                raise Warning("Something went wrong with paramfit during 'fitting parameters'!")

    # view results -------------------------------------------------------------
    prmfit.plot_energy(fitted_energy)

    # tleap (write new prmtop) -------------------------------------------------
    if args.forcefield == "gaff2":
        ff = "leaprc.gaff2"
    elif args.forcefield == "gaff":
        ff = "leaprc.gaff"
    else:
        raise IOError("Unknown force field {}".format(args.forcefield))

    leap_dir = "tleap"

    try:
        os.mkdir(leap_dir)
    except OSError:
        pass

    # append new frcmod to existing ones to create a new prmtop
    cur_frcmod = os.path.abspath("fitted_params.frcmod")
    frcmods.append(cur_frcmod)

    os.chdir(leap_dir)
    # write input file for tleap

    # write tleap input script
    if os.path.isfile("tleap.in") is False:
        with open("tleap.in", "w") as f_out:
            f_out.write("source {}\n".format(ff) +
                        "MOL = loadmol2 {}\n".format(args.mol2))

            for cfrcmod in frcmods:
                f_out.write("loadamberparams {}\n".format(cfrcmod))

            f_out.write("check MOL\n" +
                        "saveamberparm MOL {} {}\n".format("refined.prmtop", "refined.inpcrd") +
                        "quit\n")

    # write prmtop
    if os.path.isfile("refined.prmtop") is False:
        write_prmtop = sp32.call([tleap, "-f", "tleap.in"])

    cur_prmtop = os.path.abspath("refined.prmtop")

    # get back to working directory
    os.chdir(workdir)
