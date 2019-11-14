#!/usr/bin/env python

import argparse
import ag_unify_md as agum

# Argument Parsing -------------------------------------------------------------
parser = argparse.ArgumentParser(
    prog="mae2mdcrd_gau.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Convert maestro (*.mae)- to amber-trajectory (*.mdcrd)- and \
                 gaussian-input-files. Intended to be used for paramfit (Ambertools) utility.")

parser.add_argument("mae",
                    metavar="*.mae",
                    help="Maestro-file.",
                    )

parser.add_argument("-nproc",
                    default=4,
                    type=int,
                    help="Number of available CPU-Cores."
                    )

parser.add_argument("-mem",
                    default=4,
                    type=int,
                    help="Amount of available system memory in GB integers."
                    )

parser.add_argument("-job_type",
                    default="SP",
                    help="Type of job to run in gaussian " +
                         "(geometry optimization, single point energy calculation)."
                    )

parser.add_argument("-method",
                    default="MP2",
                    help="Level of theory for the gaussian job."
                    )

parser.add_argument("-basis_set",
                    default="6-311++G**",
                    help="Basis set to use by gaussian."
                    )

parser.add_argument("-geom",
                    default="connectivity",
                    help="Geometry keyword for gaussian (PrintInputOrient|connectivity)."
                    )

parser.add_argument("-charge",
                    default=0,
                    type=int,
                    help="Charge of the molecule."
                    )

parser.add_argument("-multiplicity",
                    default=1,
                    action="store",
                    type=int,
                    help="Multiplicity (2S+1) of the molecule."
                    )

parser.add_argument("-out",
                    default="DEFAULTNAME",
                    help="Name of output-files."
                    )

args = parser.parse_args()

# Read mae-file ----------------------------------------------------------------
mae = agum.Unification()
# get conformers and bonding orders from maestro's jaguar calculation
mae.read_mae(args.mae)

# frame-ids
frame_ids = list(range(len(mae.ts_coords)))

# write amber trajectory for paramfit
mae.write_mdcrd(args.out+".mdcrd", frame_ids)

# define the rest
mae.nproc        = args.nproc
mae.mem          = args.mem
mae.job_type     = args.job_type
mae.method       = args.method
mae.basis_set    = args.basis_set
mae.geom         = args.geom
mae.charge       = args.charge
mae.multiplicity = args.multiplicity

# write gaussian input-files
for frame_id in frame_ids:
    # name of checkpoint
    mae.chk = "{}_{}".format(args.out, frame_id)
    # write gaussian input file
    mae.write_gau(args.out+"_{}.gau".format(frame_id), frame_id,
                  title="Conformer {}".format(frame_id))
