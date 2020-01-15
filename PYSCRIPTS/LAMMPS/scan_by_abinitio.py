#!/usr/bin/env python
import pdb
import os
from lammps import lammps
import ag_fileio
import ag_unify_md as agum

ABINITIO_SOURCEPATH = "/hades/gadelmeier/Research.new/carbamazepine/2.ab_initio/3.scans/1.dimer_scans/0H_anti_rigid_scan/1.vertical_scan/1.gaussian/MP2_def2TZV/"
FILEENDING = "*.gau.out"

if __name__ == "__main__":
    FH = ag_fileio.FileHandler()
    FH.find_files(ABINITIO_SOURCEPATH, FILEENDING)
    FH.files_to_strings()
    MAINSYS = agum.Unification()

    for gaufile in FH.files:
        if gaufile.endswith(".gau.out"):
            MAINSYS.read_gau_log(gaufile, read_summary=True)
            pdb.set_trace()
