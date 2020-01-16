#!/usr/bin/env python
import pdb
import os
from lammps import lammps, PyLammps
import ag_fileio
import ag_unify_md as agum
import ag_lammps_sim as aglmpsim

ABINITIO_SOURCEPATH = "/hades/gadelmeier/Research.new/carbamazepine/2.ab_initio/3.scans/1.dimer_scans/0H_anti_rigid_scan/1.vertical_scan/1.gaussian/MP2_def2TZV/"
FILEENDING = ".gau.out"
PCFILE = ""
SETTINGS = ""

LMPSETTINGS = aglmpsim.LmpSim(tstart=1, tstop=1, logsteps=1, runsteps=0, pc_file=PCFILE, settings_file=SETTINGS, input_lmpdat="", output_lmplog="", )

def ab_initio2xyz(mainpath, filetype=".gau.out"):
    FH = ag_fileio.FileHandler()
    FH.find_files(ABINITIO_SOURCEPATH, FILEENDING)
    FH.files_to_strings()
    MAINSYS = agum.Unification()

    for gaufile in FH.files:
        if gaufile.endswith(".gau.out"):
            MAINSYS.read_gau_log(gaufile, read_summary=True)

    frame_ids = range(len(MAINSYS.ts_coords))
    MAINSYS.write_xyz("tmp.xyz", *frame_ids, title="DEFAULT", guess_element=False)



if __name__ == "__main__":
    pass
