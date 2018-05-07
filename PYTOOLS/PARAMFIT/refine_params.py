#!/usr/bin/env python
from __future__ import print_function
import os
import shutil as sl
import argparse
import subprocess32 as sp32
import time

__version__ = "2017-08-07"

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="refine_params.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Convert maestro (*.mae)- to amber-trajectory (*.mdcrd)- and \
                 gaussian-input-files. Used for Ambertools/paramfit utility.")

parser.add_argument("-prmtop",
                    required=True,
                    metavar="*.prmtop",
                    action="store",
                    help="Amber topology file.",
                    )

parser.add_argument("-mdcrd",
                    required=True,
                    metavar="*.prmtop",
                    action="store",
                    help="Amber coordinates file with conformations from \
                         relaxed molecule scan.",
                    )

parser.add_argument("-quantum",
                    required=True,
                    metavar="*.dat",
                    action="store",
                    help="Extracted single point energies from gaussian calculations. \
                          order must be the same as in amber coordinates file!",
                    )

parser.add_argument("-nstructures",
                    default=81,
                    type=int,
                    action="store",
                    help="Number of conformations from relaxed scan.",
                    )

parser.add_argument("-mol2",
                    required=True,
                    help="MOL2 file of the structure which is optimized. " +
                         "Must be the file that was used to parametrize the molecule in the first place!"
                    )

parser.add_argument("-frcmods",
                    nargs="*",
                    required=True,
                    help="Amber frcmod file which has parameters that were not represented properly by gaff(2)")

args = parser.parse_args()

# define environment vars ------------------------------------------------------
AMBERHOME   = os.getenv("AMBERHOME")
paramfit    = "{}/bin/paramfit".format(AMBERHOME)
antechamber = "{}/bin/antechamber".format(AMBERHOME)
tleap       = "{}/bin/tleap".format(AMBERHOME)
pwd_files   = os.listdir("./")

# 'set_params.in'; starts paramfit in interactive mode, asks for parameters to fit
# write file only if another is not present already
set_params_in  = "set_params.in"
set_params_out = "set_params.out"

if set_params_in not in pwd_files:
    with open(set_params_in, "w") as f_out:
        f_out.write(
            "RUNTYPE=SET_PARAMS\n" +
            "PARAMETER_FILE_NAME={}\n".format(set_params_out)
        )

# only write file/ run paramfit if none is present already
if set_params_out not in pwd_files:
    # set parameters using paramfit
    set_parms = sp32.call([paramfit, "-i", "set_params.in", "-p", args.prmtop])
    time.sleep(5)


# 'fit_K.in'; calculates force constant K for given parameter
fit_k_in  = "fit_K.in"
fit_k_out = "fit_K.out"

if fit_k_in not in pwd_files:
    with open(fit_k_in, "w") as f_out:
        f_out.write(
            "RUNTYPE=FIT\n"
            "PARAMETERS_TO_FIT=K_ONLY\n" +
            "NSTRUCTURES={}\n".format(args.nstructures) +
            "COORDINATE_FORMAT=TRAJECTORY\n" +
            "FUNC_TO_FIT=SUM_SQUARES_AMBER_STANDARD\n" +
            "QM_ENERGY_UNITS=HARTREE\n"
        )

if fit_k_out not in pwd_files:
    # open file for writing stdout
    with open(fit_k_out, "w") as out:
        print("***Info: Exexuting paramfit (initial K)")
        fit_K = sp32.call([paramfit, "-i", fit_k_in, "-p", args.prmtop,
                           "-c", args.mdcrd, "-q", args.quantum], stdout=out)

# extract K from fit_K.out
with open(fit_k_out, "r") as f_in:
    for line in f_in:

        if "FINAL PARAMETERS" in line:
            # skip the following two lines
            f_in.next()
            f_in.next()
            force_const_K = float(f_in.next().split()[2])
            print("***Info: Found value of K: ", force_const_K)
            break

# 'fit_params.in'; utilizes K from previous step, does the actual fit
fit_params_in  = "fit_params.in"
fit_params_out = "fit_params.out"
fit_energies   = "fit_output_energy.dat"
frcmod_fit     = "fitted_params.frcmod"

if fit_params_in not in pwd_files:
    with open(fit_params_in, "w") as f_out:
        f_out.write(
            "RUNTYPE=FIT\n" +
            "PARAMETERS_TO_FIT=LOAD\n" +
            "PARAMETER_FILE_NAME={}\n".format(set_params_out) +
            "COORDINATE_FORMAT=TRAJECTORY\n" +
            "NSTRUCTURES={}\n".format(args.nstructures) +
            "K={}\n".format(force_const_K) +
            "FUNC_TO_FIT=SUM_SQUARES_AMBER_STANDARD\n" +
            "QM_ENERGY_UNITS=HARTREE\n" +
            "ALGORITHM=BOTH\n" +
            "OPTIMIZATIONS=500\n" +
            "MAX_GENERATIONS=10000\n" +
            "GENERATIONS_TO_CONV=20\n" +
            "GENERATIONS_TO_SIMPLEX=5\n" +
            "GENERATIONS_WITHOUT_SIMPLEX=5\n" +
            "MUTATION_RATE=0.100000\n" +
            "PARENT_PERCENT=0.250000\n" +
            "SEARCH_SPACE=-1.000000\n" +
            "SORT_MDCRDS=OFF\n" +
            "WRITE_ENERGY={}\n".format(fit_energies) +
            "WRITE_FRCMOD={}\n".format(frcmod_fit)
        )

if fit_params_out not in pwd_files:
    with open(fit_params_out, "w") as out:
        print("***Info: Executing paramfit (fitting parameters)")
        fit_param = sp32.call([paramfit, "-i", "fit_params.in", "-p", args.prmtop,
                              "-c", args.mdcrd, "-q", args.quantum], stdout=out)

# write amber topology file file using tleap (and a coordinate file which is not of interest here)
with open("tleap.in", "w") as f_out:
    f_out.write("source leaprc.gaff2\n")
    f_out.write("MOL = loadmol2 {}\n".format(args.mol2))
    for cfrcmod in args.frcmods:
        f_out.write("loadamberparams {}\n".format(cfrcmod))
    f_out.write("loadamberparams {}\n".format(frcmod_fit))
    f_out.write("check MOL\n")
    f_out.write("saveamberparm MOL refined.prmtop mol2.inpcrd\n")
    f_out.write("quit\n")

write_prmtop = sp32.call([tleap, "-f", "tleap.in"])

try:
    os.mkdir("tleap")
except OSError:
    pass

sl.copyfile("refined.prmtop", "tleap/refined.prmtop")
sl.copyfile("mol2.inpcrd", "tleap/mol2.inpcrd")
sl.copyfile("tleap.in", "tleap/tleap.in")
sl.copyfile("leap.log", "tleap/leap.log")
# overwrite stuff
os.remove("refined.prmtop")
os.remove("mol2.inpcrd")
os.remove("tleap.in")
os.remove("leap.log")
