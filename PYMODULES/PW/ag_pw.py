"""
PW Module.

Read/Write PW-Input-/Output-Files.
"""

from __future__ import print_function, division
import math
import os
import re
import numpy as np
import scipy.constants as sc
# import time
import md_box as mdb
import md_stars as mds
import md_universe as mdu
# import ag_vectalg as agv

__version__ = "2018-06-22"

BOHR_ANGSTROM = sc.value("Bohr radius") / sc.angstrom
ANGSTROM_BOHR = sc.angstrom / sc.value("Bohr radius")
RYDBERG_EV = sc.value("Rydberg constant times hc in eV")


class PwStuff(mdu.Universe):
    """
    PW Stuff.

    Read and write PW input/output files.
    """

    def __init__(self):
        """
        Initialize general stuff.
        """
        # generates lists for atom-, bond-, angle-types, etc. pp.
        mdu.Universe.__init__(self)
        self.pw_entries = {}
        self.pw_entries["CONTROL"] = {}
        self.pw_entries["SYSTEM"] = {}
        self.pw_entries["ELECTRONS"] = {}
        self.pw_entries["IONS"] = {}
        self.pw_entries["CELL"] = {}
        self.pw_entries["K_POINTS"] = {}
        self.pw_other_info = {}
        self.pw_other_info["ENERGIES"] = []
        self.pw_other_info["DENSITIES"] = []
        self.pw_other_info["VOLUMES"] = []

    def read_pwin(self, pwin):
        """
        Read an input file for pw.x.

        Cell vector alat (celldm(1)) is converted to angstrom when read.
        """
        with open(pwin) as opened_pwin:
            line = opened_pwin.readline()
            while line != '':

                if line.startswith("&CONTROL"):
                    while line != '':
                        line = opened_pwin.readline()

                        # skip empty lines, end if end of block ('/') is reached
                        if line.startswith("\n"):
                            continue
                        elif line.startswith("/"):
                            break
                        else:
                            pass

                        # split line by '=' and ' '
                        split_line = re.split(r'\s*=\s*', line)
                        split_line[1] = split_line[1].strip()

                        # convert str float/int if possible
                        if re.match(r"^\d+?\.\d+?$", split_line[1]):
                            split_line[1] = float(split_line[1])
                        else:
                            try:
                                split_line[1] = int(split_line[1])
                            except ValueError:
                                pass

                        self.pw_entries["CONTROL"][split_line[0].strip()] = split_line[1]
                elif line.startswith("&SYSTEM"):
                    while line != '':
                        line = opened_pwin.readline()

                        if line.startswith("\n"):
                            continue
                        elif line.startswith("/"):
                            break
                        else:
                            pass

                        split_line = re.split(r'\s*=\s*', line)
                        split_line[1] = split_line[1].strip()

                        # convert str float/int if possible
                        if re.match(r"^\d+?\.\d+?$", split_line[1]):
                            split_line[1] = float(split_line[1])
                        else:
                            try:
                                split_line[1] = int(split_line[1])
                            except ValueError:
                                pass

                        self.pw_entries["SYSTEM"][split_line[0].strip()] = split_line[1]
                elif line.startswith("&ELECTRONS"):
                    while line != '':
                        line = opened_pwin.readline()

                        if line.startswith("\n"):
                            continue
                        elif line.startswith("/"):
                            break
                        else:
                            pass

                        split_line = re.split(r'\s*=\s*', line)
                        split_line[1] = split_line[1].strip()

                        # convert str float/int if possible
                        if re.match(r"^\d+?\.\d+?$", split_line[1]):
                            split_line[1] = float(split_line[1])
                        else:
                            try:
                                split_line[1] = int(split_line[1])
                            except ValueError:
                                pass

                        self.pw_entries["ELECTRONS"][split_line[0].strip()] = split_line[1]
                elif line.startswith("&IONS"):
                    while line != '':
                        line = opened_pwin.readline()

                        if line.startswith("\n"):
                            continue
                        elif line.startswith("/"):
                            break
                        else:
                            pass

                        split_line = re.split(r'\s*=\s*', line)
                        split_line[1] = split_line[1].strip()

                        # convert str float/int if possible
                        if re.match(r"^\d+?\.\d+?$", split_line[1]):
                            split_line[1] = float(split_line[1])
                        else:
                            try:
                                split_line[1] = int(split_line[1])
                            except ValueError:
                                pass

                        self.pw_entries["IONS"][split_line[0].strip()] = split_line[1]
                elif line.startswith("&CELL"):
                    while line != '':
                        line = opened_pwin.readline()

                        if line.startswith("\n"):
                            continue
                        elif line.startswith("/"):
                            break
                        else:
                            pass

                        split_line = re.split(r'\s*=\s*', line)
                        split_line[1] = split_line[1].strip()

                        # convert str float/int if possible
                        if re.match(r"^\d+?\.\d+?$", split_line[1]):
                            split_line[1] = float(split_line[1])
                        else:
                            try:
                                split_line[1] = int(split_line[1])
                            except ValueError:
                                pass

                        self.pw_entries["CELL"][split_line[0].strip()] = split_line[1]
                elif line.startswith("ATOMIC_SPECIES"):
                    # GET ATOM TYPES
                    atmcnt = 0

                    while line != '':
                        line = opened_pwin.readline()

                        # end of block
                        if line == "\n":
                            break

                        split_line = line.split()
                        # add current atom type to atom types
                        self.atm_types[atmcnt] = mds.Atom(sitnam=split_line[0],
                                                          weigh=float(split_line[1]),
                                                          pseudopotential=split_line[2])
                        atmcnt += 1
                elif line.startswith("ATOMIC_POSITIONS"):
                    # GET ATOM COORDINATES
                    self.ts_coords.append([])

                    while line != '':
                        line = opened_pwin.readline()

                        # end of block
                        if line == "\n":
                            break

                        split_line = line.split()
                        self.atoms.append(mds.Atom(sitnam=split_line[0]))
                        # add coordinates from current atom to the current frame
                        self.ts_coords[-1].append(np.array([float(i) for i in split_line[1:]]))

                elif line.startswith("K_POINTS"):
                    kpoints_line = line.split()
                    self.pw_entries["K_POINTS"]["option"] = kpoints_line[1]

                    # read upcoming lines
                    while line != '':
                        line = opened_pwin.readline()

                        if line.startswith("\n"):
                            break

                        split_line = line.split()

                        if kpoints_line[1].strip("{}()") == "automatic":
                            self.pw_entries["K_POINTS"]["k_point_grid"] = [int(i) for i in split_line]
                            break
                        elif kpoints_line[1] == "gamma":
                            self.pw_entries["K_POINTS"]["k_point_grid"] = []
                            break
                        else:
                            pass

                elif line.startswith("CELL_PARAMETERS"):
                    box_unit = line.split()[1].strip("{}")

                    # get box vectors
                    cbox = mdb.Box(boxtype="cartesian", unit=box_unit)

                    for line_cntr in xrange(3):
                        line = [float(i) for i in opened_pwin.readline().split()]

                        # allot each vector to the box
                        if line_cntr == 0:
                            cbox.crt_a = line  # vector a
                        elif line_cntr == 1:
                            cbox.crt_b = line  # vector b
                        else:
                            cbox.crt_c = line  # vector c

                    self.ts_boxes.append(cbox)
                else:
                    pass

                line = opened_pwin.readline()

        # convert cell to lattice cell and append it to the current cells
        if ("A" in self.pw_entries["SYSTEM"] and
            "B" in self.pw_entries["SYSTEM"] and
            "C" in self.pw_entries["SYSTEM"] and
            "cosAB" in self.pw_entries["SYSTEM"] and
            "cosAC" in self.pw_entries["SYSTEM"] and
            "cosBC" in self.pw_entries["SYSTEM"]):
            #
            cbox = mdb.Box(
                ltc_a=float(self.pw_entries["SYSTEM"]["A"]),
                ltc_b=float(self.pw_entries["SYSTEM"]["B"]),
                ltc_c=float(self.pw_entries["SYSTEM"]["C"]),
                ltc_alpha=math.acos(float(self.pw_entries["SYSTEM"]["cosBC"])),
                ltc_beta=math.acos(float(self.pw_entries["SYSTEM"]["cosAC"])),
                ltc_gamma=math.acos(float(self.pw_entries["SYSTEM"]["cosAB"])),
                boxtype="lattice",
                unit="angstrom")

            # delete surplus box entries
            del self.pw_entries["SYSTEM"]["A"]
            del self.pw_entries["SYSTEM"]["B"]
            del self.pw_entries["SYSTEM"]["C"]
            del self.pw_entries["SYSTEM"]["cosBC"]
            del self.pw_entries["SYSTEM"]["cosAC"]
            del self.pw_entries["SYSTEM"]["cosAB"]

            # add celldm(1) and convert it to bohr
            #self.pw_entries["SYSTEM"]["celldm(1)"] = cbox.ltc_a * ANGSTROM_BOHR
            self.pw_entries["SYSTEM"]["ibrav"] = 0

            # convert lattice box to cartesian
            cbox.box_lat2cart()
            self.ts_boxes.append(cbox)

        # convert celldm (= alat) to box vector a with angstrom
        #try:
            #self.ts_boxes[-1].ltc_a = float(self.pw_entries["SYSTEM"]["celldm(1)"]*BOHR_ANGSTROM)
        #except KeyError:
            #pass

        # final check
        if len(self.atoms) != self.pw_entries["SYSTEM"]["nat"]:
            print("***Warning: Number of atoms and SYSTEM entry 'nat' differ!")
            #time.sleep(5)

        if os.path.isdir(self.pw_entries["CONTROL"]["pseudo_dir"].strip("'")) is False:
            print("***Warning: Folder for Pseudopotentials does not exist!")
            #time.sleep(5)

    def read_pwout(self, pwout):
        """
        CAVEAT: UNDER CONSTRUCTION! Read the output of pw.x.

        Currently this only reads the coordinates and cell vectors.
        Cell vector alat (celldm(1)) is converted to angstrom when read.
        """
        #print(pwout)
        with open(pwout) as opened_pwout:
            line = opened_pwout.readline()
            while line != '':

                if line.startswith("CELL_PARAMETERS"):
                    # get alat
                    split_line = line.split()

                    # get box vectors
                    cbox = mdb.Box(boxtype="cartesian")

                    if "alat" in line:
                        #self.pw_entries["SYSTEM"]["celldm(1)"] = float(split_line[2].strip(")"))
                        cbox.unit = "alat"
                    elif "bohr" in line:
                        cbox.unit = "bohr"
                    elif "angstrom" in line:
                        cbox.unit = "angstrom"
                    else:
                        raise Warning("Keyword for box unknown and not implemented.")

                    cbox.crt_a = [float(i)
                                  for i in opened_pwout.readline().split()]
                    cbox.crt_b = [float(i)
                                  for i in opened_pwout.readline().split()]
                    cbox.crt_c = [float(i)
                                  for i in opened_pwout.readline().split()]

                    if cbox.unit == "alat":
                        cbox.alat2angstrom(float(split_line[2].strip(")")))

                    self.ts_boxes.append(cbox)

                elif line.startswith("ATOMIC_POSITIONS"):
                    self.atoms = []  # overwrite existing atoms
                    # prepare container for coordinates to come
                    self.ts_coords.append([])

                    # read the coordinates
                    while line != '':
                        line = opened_pwout.readline()

                        # stop reading when EOF is reached
                        if line.startswith("\n") or line.startswith("End final coordinates"):
                            break

                        split_line = line.split()

                        cur_atm = mds.Atom(sitnam=split_line[0])
                        self.atoms.append(cur_atm)
                        cur_atm_coords = np.array([float(i) for i in split_line[1:]])
                        self.ts_coords[-1].append(cur_atm_coords)
                elif line.startswith("!    total energy"):
                    line = line.split()
                    energy = float(line[-2]) * RYDBERG_EV
                    self.pw_other_info["ENERGIES"].append(energy)
                elif line.startswith("     density ="):
                    line = line.split()
                    density = float(line[-2])
                    self.pw_other_info["DENSITIES"].append(density)
                elif line.startswith("     new unit-cell volume"):
                    line = line.split()
                    volume = float(line[-3])
                    self.pw_other_info["VOLUMES"].append(volume)
                else:
                    pass

                line = opened_pwout.readline()

    def _write_section(self, opened_file_instance, frame_id, keyword):
        """
        Help writing an input section.

        opened_file_instance    file; file to write section to
        frame_id                int; id of coordinates to write
        keyword                 str; pw_entry keyword
        """
        # write entries
        opened_file_instance.write("&{}\n".format(keyword))

        for setting, value in self.pw_entries[keyword].iteritems():

            # calculate celldm(1), which derives from lattice vector a
            if setting == "celldm(1)":
                opened_file_instance.write("    {0:<24s}= {1}\n".format(
                    setting, self.ts_boxes[frame_id].ltc_a * ANGSTROM_BOHR))
            else:
                opened_file_instance.write("    {0:<24s}= {1}\n".format(
                    setting, value))

        opened_file_instance.write("/\n\n")

    def write_pwin(self, frame_id, filename, verbosity="'high'"):
        """
        Write an input file for pw.x. alat is converted to angstrom when read.

        frame_id                int; id of coordinates to write
        filename                str; pw-output file to write to
        verbosity               str; high | low
        """
        self.pw_entries["CONTROL"]["verbosity"] = verbosity
        #self.pw_entries["SYSTEM"]["A"] = agv.get_mag(self.ts_boxes[frame_id].crt_a)

        with open(filename, "w") as opened_filename:
            # since it seems that a certain order must be maintained,
            # section by section is written

            # &CONTROL Section
            self._write_section(opened_filename, frame_id, "CONTROL")
            self._write_section(opened_filename, frame_id, "SYSTEM")
            self._write_section(opened_filename, frame_id, "ELECTRONS")
            self._write_section(opened_filename, frame_id, "IONS")
            self._write_section(opened_filename, frame_id, "CELL")

            # Pseudopotentials (atom types)
            opened_filename.write("ATOMIC_SPECIES\n")
            for _, atom_type in self.atm_types.iteritems():
                opened_filename.write("{0} {1:>6} {2:}\n".format(
                    atom_type.sitnam, atom_type.weigh,
                    atom_type.pseudopotential))
            opened_filename.write("\n")

            # Atom Coordinates
            opened_filename.write("ATOMIC_POSITIONS {angstrom}\n")
            coordinate_string = "{0:<5s} {c[0]:> 18.9f}{c[1]:> 18.9f}{c[2]:> 18.9f}\n"
            for atom, atom_coordinates in zip(self.atoms,
                                              self.ts_coords[frame_id]):
                opened_filename.write(coordinate_string.format(
                    atom.sitnam, c=atom_coordinates))

            opened_filename.write("\n")

            # K Points
            opened_filename.write("K_POINTS {}\n".format(
                self.pw_entries["K_POINTS"]["option"]))

            # write grid or skip entry if gamma point is set
            if self.pw_entries["K_POINTS"]["option"].strip("{}()") == "automatic":
                for kpoint in self.pw_entries["K_POINTS"]["k_point_grid"]:
                    opened_filename.write("{} ".format(kpoint))
                opened_filename.write("\n")

            opened_filename.write("\n")

            # Cell
            # write section only if ibrav is 0
            try:
                if self.pw_entries["SYSTEM"]["ibrav"] == 0:
                    cell_string = "{v[0]:> 24.9f}{v[1]:> 18.9f}{v[2]:> 18.9f}\n"
                    opened_filename.write("CELL_PARAMETERS {{{0}}}\n".format(
                        self.ts_boxes[frame_id].unit))
                    opened_filename.write(
                        cell_string.format(v=self.ts_boxes[frame_id].crt_a))
                    opened_filename.write(
                        cell_string.format(v=self.ts_boxes[frame_id].crt_b))
                    opened_filename.write(
                        cell_string.format(v=self.ts_boxes[frame_id].crt_c))
            except KeyError:
                pass


def read_pw_out(pw_out):
    """Read a pw output file."""
    pw_sys = PwStuff()
    pw_sys.read_pwout(pw_out)
    return pw_sys
