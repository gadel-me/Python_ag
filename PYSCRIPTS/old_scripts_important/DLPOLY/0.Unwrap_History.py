#!/home/angad/Programs/anaconda/bin/python -i

from __future__ import print_function, division
import os
import tqdm
import numpy as np
import sys
sys.path.append("../myPYMODULES/DLPOPLY_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/OTHER_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/MATH_MODULES/")
import DL_HISTORY_class8_new_approaches_8_2 as dlhc
#import my_lin_alg_compendium as mlc
import argparse

# Argument parsing stuff
parser = argparse.ArgumentParser(prog="DL_POLY CONFIG/HISTORY rotator and unwrapper",
                                 description="Rotate and/or unwrap a DL_POLY CONFIG and/or DL_POLY HISTORY-File")

parser.add_argument("-c", "--dlpolyconfig",
                    help="DL_POLY CONFIG or REVCON-File with box vectors",
                    action="store")

parser.add_argument("-cout", "--dlpolyconfig_out",
                    help="Filename of modified CONFIG (default='DEFAULTNAME.dlpolyconfig'",
                    action="store")

parser.add_argument("-hout", "--dlpolyhist_out",
                    help="Filename of modified HISTORY (default='DEFAULTNAME.dlpolyhist'",
                    action="store")

parser.add_argument("-h", "--dlpolyhist",
                    help="DL_POLY HISTORY-File with box vectors for each time step",
                    action="store")

parser.add_argument("-f", "--dlpolyfld",
                    required=True,
                    help="DL_POLY FIELD-File belonging to the HISTORY/CONFIG/REVCON-File",
                    action="store")

parser.add_argument("--rotate",
                    help="Rotate the whole box so it can be understood by VMD and LAMMPS.\nBox-Vector a: [x, 0, 0]\nBox-Vector b:[x, y, 0]\nBox-Vector c:[x, y, +z]",
                    action="store_true")

parser.add_argument("--unwrap",
                    help="Unwrap the whole box.\nAtoms will be outside of periodic boundaries.\n***Note: It is possible that Molecules get shifted unreasonably.)",
                    action="store_true")

parser.add_argument("--step_start",
                    help="Start processing the HISTORY-File from this timestep",
                    type=int,
                    default=None)

parser.add_argument("--step_stop",
                    help="Stop processing the HISTORY-File beyond this timestep",
                    type=int,
                    default=None)

parser.add_argument("--version",
                    help="Version information",
                    action="store_true",
                    version=0.90)

args = parser.parse_args()

script, historyfile, configfile, fieldfile, step_start, step_stop = sys.argv

#TODO renew description
"""
The aim of this program is to calculate the shifts of the mass centers and the
intra-molecular angles (relative to the "original" vectors).
Original vectors descend from the crystallographic data (CIF-Files), which
were converted to the DL_POLY CONFIG-format. Using these files, the center of
masses and the intra-molecular coordinate systems are calculated for each
molecule. The base vectors of these intra molecular coordinate systems are
the reference for the angles chi (x-vector (cif) vs. x-vector (MD)), psi
(y-vector (cif) vs. y-vector (MD) and omega (z-vector (cif) vs. z-vector (MD)).

> moltype:  molecule name in the FIELD (e.g. TIP3P)
> COM:      center of mass
> SUM_COMs: summation of all center of masses
> AXIS_X:   intra-molecular x-axis
> AXIS_Y:   intra-molecular y-axis
> AXIS_Z:   intra-molecular z-axis
> CHI:      angle between intra-x-axis and base x
> PSI:      angle between intra-y-axis and base y
> OMEGA:    angle between intra-z-axis and base z
"""
# DELETE THE REMNANTS OF PREVIOUS ATTEMPTS /////////////////////////////////////
if os.path.isfile("HISTORY.dlpolyhist"):
    os.remove("HISTORY.dlpolyhist")

# PROCESS ONLY STEPS STARTING FROM 'step_start' ////////////////////////////////
# 'step_start' = step/'traj-value' (from CONTROL-File)
step_start = int(step_start)

# FIELD-INSTANCE ///////////////////////////////////////////////////////////////
finfo = dlhc.FieldFile(fieldfile)
finfo.read_field()

# CONFIG = REFERENCE ///////////////////////////////////////////////////////////
ref_conf = dlhc.ConfigFile(finfo, configfile)
ref_conf.read_config()
# COM-SHIFT => PUSH UNWRAP TO CERTAIN DIRECTION
ref_conf.config.unwrap_box(
    finfo,
    mass_center_shift=np.array([10, 10, 10])
)
ref_conf.config.write_config()
# INTRAMOLECULAR COORDINATE SYSTEM = ROTATIONAL DEGREES OF FREEDOM
##ref_conf.config.molecule_types[0].get_all_mcs(25, 26, 13, 11)
# CENTER OF MASSES FOR THE CURRENT MOLECULE-POSITIONS
##ref_conf.config.molecule_types[0].get_all_coms()
# GIVE REFERENCE A NEW TAG FOR UNWRAPPING THE WHOLE HISTORY (see below)
cur_conf_ref = ref_conf.config
# convert all config-mcs types ('object' to 'float64')
for i in ref_conf.config.molecule_types[0].molecules:
    i.coord_sys = [j.astype("float64") for j in i.coord_sys]

# STEP 1 ///////////////////////////////////////////////////////////////////////
# unfold each box from each frame according to the reference
hist = dlhc.HistoryFile(finfo, historyfile)
hist.read_nsteps(ts_stop=None, debug=False)
# better than re-reading file after every instance of 'HistoryFile' all over again
opened_history = open(historyfile)
# get all steps that are in the history (performance reasons, see below)
process_keys = hist.nstep_keys[step_start:]
unwrapped_history = "DEFAULT.dlpolyhist"
# stuff for making durance visible
pbar = tqdm.tqdm(process_keys)

# large history files must be read this way because the memory-overhead is
# be too massive
for istep, tqdmstep in zip(process_keys, pbar):
    # creating new instances of the class "HistoryFile" makes it possible
    # to fit large files in memory
    cur_hstep = dlhc.HistoryFile(finfo, opened_history)
    cur_hstep.read_history(ts_start=istep, ts_stop=istep)
    # unwrap each box according to the reference
    cur_hstep.timesteps[istep].config.unwrap_box(finfo)
    # "max_discr" is the maximal allowed distance between the two
    # center of masses that are compared currently; makes the comparison run
    # faster because not every molecule has to be tested against each other
    cur_hstep.timesteps[istep].config.shift_all_mols_to_ref(
        cur_conf_ref,
        max_discr=10)
    # rotate box so vmd can understand it
    cur_hstep.timesteps[istep].config.rotate_box_for_vmd()
    # every step is appended
    cur_hstep.write_history(history_out=unwrapped_history)

    # correct current position of file:
    # method "read_history" does not stop until having read the following
    # "TIMESTEP"-line -> file must be rewound to a point before that line

    # get current position
    cur_hist_pos = opened_history.tell()
    # cur_hist_pos must be large enough cause file-read is buffered
    cur_hist_pos -= 100000
    # rewind by cur_hist_pos (to ensure to get before the "TIMESTEP"-line)
    opened_history.seek(cur_hist_pos)
    # redefine the reference (to be the current config)
    cur_conf_ref = cur_hstep.timesteps[istep].config
    # stuff for making durance visible
    pbar.set_description("Processing %s" % istep)
    pbar.update(1)

opened_history.close()
