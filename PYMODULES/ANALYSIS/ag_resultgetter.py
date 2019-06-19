"""
This module is intended for use with jupyter notebook 'Notebook - 1'.
It provides functions to find gaussian output, pw-output and lammps-log files,
extract their geometries with the according energies and measures the distances
between the scanned molecules by given indices.

Examples
--------

Gaussian and Quantum Espresso output files (Measure distance between two center of geometries and read the according energy)
------------------------------------------
results = ResultGetter("path/to/resultfolder", "*.gau.out")
results.process_pw_gau(
    [14, 15, 16, 18, 20, 22],
    [44, 45, 46, 48, 50, 52],
    "gau_out")

Lammps log and dcd files
------------------------
results = ResultGetter()
results.process_lmp(
    "lmp/log/file",
    "dcd/file",
    [14, 15, 16, 18, 20, 22],
    [44, 45, 46, 48, 50, 52])

"""

import pdb
import os
import numpy as np
import progressbar
import ag_fileio
import ag_unify_md as agum
import ag_unify_log as agul


class ResultGetter(ag_fileio.FileHandler):
    """
    """

    def __init__(self, directory=None, name_pattern=None):
        """
        Immediately search for files.
        """
        ag_fileio.FileHandler.__init__(self)

        # find files according to name_pattern
        if directory is not None and name_pattern is not None:
            self.find_files(directory, name_pattern)
            #self.files = [os.path.abspath(i) for i in self.files]

        # list of results
        self.energy_unit = "eV"
        self.results = []
        self.normed_results = []

    def _norm_results(self):
        # sort results by distance (first value in tuple)
        sorted_results = sorted(self.results, key=lambda x: float(x[0]))
        # subtract the energy from the furthest distance between both dimers
        substr_energy = sorted_results[-1][-1]
        # normalize all energies against above named energy
        self.normed_results = [[i[0], (i[1] - substr_energy)] for i in sorted_results]

    @staticmethod
    def _get_distance(mdsys, idxs1, idxs2, frame_id=-1):
        """
        Calculate the distance between two atoms or between the cogs of two atom-clusters.
        """
        if len(idxs1) == 1:
            cog1 = mdsys.ts_coords[frame_id][idxs1[0]]
        else:
            cog1 = mdsys.get_cog(frame_id, *idxs1)

        if len(idxs2) == 1:
            cog2 = mdsys.ts_coords[frame_id][idxs2[0]]
        else:
            cog2 = mdsys.get_cog(frame_id, *idxs2)

        distance = np.linalg.norm(cog1 - cog2)
        return distance

    def _dist_and_energy_single_geom(self, filename, idxs1, idxs2, filetype=None):
        """
        Read the distance and the energy from a single output file.

        This method only works with the following setup:
            >   Each calculation has only ONE geometry to optimize and/or to
                calculate the energy from
            >   if more geometries are scanned (in a single output file), only
                the last one will be evaluated, since only the last frame is
                processed

        Parameters
        ----------
        filename : str
            name of the file to process

        idxs1 : tuple or list of ints
            atom indices to form the first cog

        idxs2 : tuple or list of ints
            atom indices to form the second cog

        filetype : bool
            type of file to read ('.pwscf_out' or '.gau_out')

        """
        dimer_sys = agum.Unification()
        # convert filename to str since PosixPath does not have an 'endswith' method
        filename = str(filename)

        if filename.endswith(".pwscf_out") or filetype == "pwscf_out":
            dimer_sys.read_pwout(filename, read_crystal_sections=True)
            energy = dimer_sys.pw_other_info["ENERGIES"][-1]
        elif filename.endswith(".gau_out") or filetype == "gau_out":
            dimer_sys.read_gau_log(filename, read_summary=True)
            energy = dimer_sys.gaussian_other_info["Counterpoise corrected energy"]
        else:
            raise IOError("File of unknown type {}".format(filetype))

        distance = self._get_distance(dimer_sys, idxs1, idxs2)
        return (distance, energy)

    def process_lmp(self, lmplog, dcd, idxs1, idxs2, frame=-1):
        """
        Extract the distance and energy from lammps output done by using 'displace_atoms' command in lammps.

        This method is reading ONE output file for the energy and for the coordinates
        which has several geometries that were optimized or single point calculations
        were carried out on.
        """
        log_data = agul.LogUnification()
        log_data.read_lmplog(lmplog)
        dimer_sys = agum.Unification()
        dimer_sys.import_dcd(dcd)
        dimer_sys.read_frames()
        nframes = xrange(len(dimer_sys.ts_coords))

        for frame_id in nframes:
            distance = self._get_distance(dimer_sys, idxs1, idxs2, frame_id)
            energy = log_data.data[frame_id]["PotEng"][frame]
            self.results.append((distance, energy))

        self._norm_results()

    def process_pw_gau(self, idxs1, idxs2, filetype=None):
        """
        Get the energies and distances from a list of pw or gaussian output files.
        """
        widgets = [
            'Processed: ',
            progressbar.Counter(),
            ' files (', progressbar.Timer(),
            ')', progressbar.Bar()
        ]
        pbar = progressbar.ProgressBar(widgets=widgets)

        for cfile in pbar(self.files):
            cresult = self._dist_and_energy_single_geom(cfile, idxs1, idxs2, filetype=filetype)
            self.results.append(cresult)

        #pbar.finish()
        self._norm_results()

    def ev_to_kcal_mol(self):
        ev_kcal = 23.061

        # convert eV to kcal/mol
        if self.energy_unit == "eV":
            self.normed_results = [[i[0], i[1] * ev_kcal] for i in self.normed_results]
            self.results = [[i[0], i[1] * ev_kcal] for i in self.results]
            self.energy_unit = "kcal/mol"
