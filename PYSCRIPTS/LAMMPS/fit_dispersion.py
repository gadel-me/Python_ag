#!/usr/bin/env python

import os
import shutil as sl
import numpy as np
import argparse
import ag_unify_md as agum
import Transformations as cgt

"""
Convert a gaussian output-/log-file (very last entry) to a xyz-file.
Do a coordinate scan between a pair of two molecules and push them away from
each other.
"""


def frange(start, stop, step):
    x = start
    while x < stop:
        yield x
        x += step


# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument(
    "-gau_out", default=None, metavar="*.lmpdat", help="Gaussian log/output-file."
)

parser.add_argument(
    "-gau_in",
    default=None,
    metavar="*.lmpdat",
    help="Gaussian input-file, suitable for having additional data.",
)

parser.add_argument(
    "-top",
    default=None,
    metavar="*.lmpdat",
    help="File with force field information and molecular topology, "
    + "e.g. lammps data file.",
)

parser.add_argument(
    "-restart", action="store_true",
)

parser.add_argument(
    "-gau_argument_line",
    default=None,
    metavar="*.lmpdat",
    help="Gaussian log/output-file.",
)

parser.add_argument(
    "-o", default=None, metavar="*.lmpdat", help="Lammps' data-file of the main system."
)

args = parser.parse_args()
# get coordinates from last gaussian log file
dimeric_sys = agum.Unification()

dimeric_sys.read_lmpdat(args.top)

# read a gaussian input file for basic information
if args.gau_in is not None:
    dimeric_sys.read_gau(args.gau_in, overwrite=True)

if args.gau_out is not None:
    dimeric_sys.read_gau_log(args.gau_out, overwrite=False)

# define job settings
if args.gau_argument_line is None:
    dimeric_sys.job_settings = "#P SP MP2(FullDirect, FC)/aug-cc-pVTZ geom=connectivity"
else:
    dimeric_sys.job_settings = args.gau_argument_line

if args.restart:
    dimeric_sys.charge = ""
    dimeric_sys.multiplicity = ""
else:
    dimeric_sys.charge = 0
    dimeric_sys.multiplicity = 1

dimeric_sys.change_indices()

if args.restart is True:
    dimeric_sys.ts_coords = []

# define hbond to shift
if args.restart is False:
    vt_shift = dimeric_sys.ts_coords[-1][25] - dimeric_sys.ts_coords[-1][55]
    unit_vt_shift = vt_shift / np.linalg.norm(vt_shift)

# delete existing pair types, write pair types to a separate file (needed when
# using dreiding h-bond potential)
# dimeric_sys.pair_types = []
# dimeric_sys.mix_pair_types(mode="ij", to_file="foo.pair")
dimeric_sys.gaussian_charges = [0]
dimeric_sys.gaussian_multiplicities = [1]
# dimeric_sys.gaussian_charges = [0, 0, 0]
# dimeric_sys.gaussian_multiplicities = [1, 1, 1]

# for i in xrange(1):
for i in frange(-2, 80, 0.5):

    output_name = "{}_{}".format(args.o, i)

    if args.restart:
        subdir = "CBZ_Dimer_anti_2H_flexible_scan_B97D3_cc-pVDZ_{}".format(i)
        dimeric_sys.oldchk = "CBZ_Dimer_anti_2H_flexible_scan_B97D3_cc-pVDZ_{}.chk".format(
            i
        )
    else:
        subdir = output_name

    os.mkdir(subdir)
    # calculate the current shift, multiply by -1 to make the shift towards
    # each other (x-coordinate is negative and we are moving according
    # the x-axis)
    if args.restart is False:
        current_shift = unit_vt_shift * (i * 0.1 * -1)
        cur_Tm = cgt.translation_matrix(current_shift)
        # translate the coordinates but let the coordinates from frame 0 stay
        # the same
        coords = dimeric_sys.mm_atm_coords(
            0, cur_Tm, True, *list(range(30, len(dimeric_sys.atoms)))
        )
        dimeric_sys.ts_coords.append(coords)
    dimeric_sys.chk = "{}.chk".format(output_name)
    dimeric_sys.write_gau(
        "{}.gau".format(output_name),
        -1,
        False,
        title="CBZ anti-dimer h-bonds scan - interatomic distance: {}".format(i),
    )

    # if args.restart is False:
    #    dimeric_sys.write_lmpdat("{}.lmpdat".format(output_name), -1,
    #                             title="CBZ anti-dimer h-bonds scan - interatomic distance: {}".format(i), cgcmm=True)

    sl.move("{}.gau".format(output_name), subdir)

    with open("start_gau.sh", "w") as sbatch_file:
        sbatch_file.write("#!/bin/bash\n")
        sbatch_file.write("#SBATCH --exclusive\n")
        sbatch_file.write("#SBATCH --job-name={}.gau\n".format(output_name))
        sbatch_file.write("g09 {}.gau".format(output_name))
    sl.move("start_gau.sh", subdir)
