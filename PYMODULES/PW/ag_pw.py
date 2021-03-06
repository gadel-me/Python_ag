"""
PW Module.

Read/Write PW-Input-/Output-Files.
When reading the output file, a common function might be used with decorators
for reading pwin-coords entry and pwout-coords entry
"""


import math
import os
import re
import pdb
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
        # container to supply atoms from "ATOMIC_POSITIONS" with info from "ATOMIC_SPECIES"
        # which may be compared to lammps "Masses" and "Atoms" input
        #TODO write frozen atoms
        atm_types_ptrs = {}

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
                        # "C": 0, "H": 1, and so on
                        atm_types_ptrs[split_line[0]] = atmcnt
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
                        # omit further information if dictionary is empty (should not happen)
                        cur_atm_sitnam = split_line[0]

                        if atm_types_ptrs == {}:
                            self.atoms.append(mds.Atom(sitnam=cur_atm_sitnam))
                        else:
                            cur_atm_key = atm_types_ptrs[cur_atm_sitnam]
                            #cur_atm_key = self.atm_types[cur_atm_key_idx]
                            catom = mds.Atom(sitnam=cur_atm_sitnam, atm_key=cur_atm_key)

                            # read frozen info if given
                            if len(split_line) == 7:
                                catom.ifrz_x = int(split_line[-3])
                                catom.ifrz_y = int(split_line[-2])
                                catom.ifrz_z = int(split_line[-1])

                            self.atoms.append(catom)

                        # add coordinates from current atom to the current frame
                        self.ts_coords[-1].append(np.array([float(i) for i in split_line[1:4]]))

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
                    #TODO: calculate from Bohr to Angstrom in place?
                    box_unit = line.split()[1].strip("{}")

                    # get box vectors
                    cbox = mdb.Box(boxtype="cartesian", unit=box_unit)

                    for line_cntr in range(3):
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

    def read_pwout(self, pwout, read_crystal_sections=False, save_all_scf_steps=True):
        """
        CAVEAT: UNDER CONSTRUCTION! Read the output of pw.x.

        Currently this only reads the coordinates and cell vectors.
        Cell vector alat (celldm(1)) is converted to angstrom when read.
        #TODO READ FROZEN ATOMS
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

                elif line.startswith("     atomic species"):
                    line = opened_pwout.readline()
                    atm_types_ptrs = {}
                    atm_type_cntr = 0

                    while line != '\n':
                        #print(repr(line))
                        atm_types_ptrs[line.split()[0]] = atm_type_cntr
                        atm_type_cntr += 1
                        line = opened_pwout.readline()

                elif line.startswith("ATOMIC_POSITIONS"):
                    #TODO: need a smarter way to do this!
                    # overwrite existing atoms
                    self.atoms = []
                    # prepare container for coordinates to come
                    self.ts_coords.append([])
                    atm_cntr = 0

                    # read the coordinates
                    while line != '':
                        line = opened_pwout.readline()

                        # stop reading when EOF is reached
                        if line.startswith("\n") or line.startswith("End final coordinates"):
                            break

                        split_line = line.split()

                        cur_atm = mds.Atom(
                            sitnam=split_line[0],
                            atm_id=atm_cntr,
                            atm_key=atm_types_ptrs[split_line[0]])
                        self.atoms.append(cur_atm)
                        cur_atm_coords = np.array([float(i) for i in split_line[1:]])
                        self.ts_coords[-1].append(cur_atm_coords)
                        atm_cntr += 1

                elif line.startswith("!    total energy"):
                    split_line = line.split()
                    energy = float(split_line[-2]) * RYDBERG_EV
                    self.pw_other_info["ENERGIES"].append(energy)
                elif line.startswith("     density ="):
                    split_line = line.split()
                    density = float(split_line[-2])
                    self.pw_other_info["DENSITIES"].append(density)
                elif line.startswith("     new unit-cell volume"):
                    split_line = line.split()
                    volume = float(split_line[-3])
                    self.pw_other_info["VOLUMES"].append(volume)
                else:
                    pass

                # get alat value
                if read_crystal_sections is True:

                    if line.startswith("     lattice parameter (alat)"):
                        # not sure why alat was converted here before, seems to be wrong
                        # when the whole box is converted later anyway
                        alat = float(line.split()[-2])  #* BOHR_ANGSTROM
                        #print(alat)

                    # get box with box vectors
                    elif line.startswith("     crystal axes: (cart. coord. in units of alat)"):
                        # get box vectors
                        cbox = mdb.Box(boxtype="cartesian")
                        cbox.unit = "alat"

                        cbox.crt_a = [float(i) for i in opened_pwout.readline().split()[3:6]]
                        cbox.crt_b = [float(i) for i in opened_pwout.readline().split()[3:6]]
                        cbox.crt_c = [float(i) for i in opened_pwout.readline().split()[3:6]]
                        cbox.alat2angstrom(alat)
                        self.ts_boxes.append(cbox)

                    elif line.startswith("   Cartesian axes"):
                        # overwrite existing atoms
                        self.atoms = []
                        # prepare container for coordinates to come
                        self.ts_coords.append([])
                        atm_cntr = 0

                        # skip next two lines
                        opened_pwout.readline()
                        opened_pwout.readline()

                        # read the coordinates
                        while line != '':
                            line = opened_pwout.readline()

                            # stop reading when end of current entry is reached
                            if line.startswith("\n"):
                                break

                            split_line = line.split()

                            cur_atm = mds.Atom(
                                sitnam=split_line[1],
                                atm_id=atm_cntr,
                                atm_key=atm_types_ptrs[split_line[1]])
                            self.atoms.append(cur_atm)
                            cur_atm_coords = np.array([float(i) * alat for i in split_line[6:9]])
                            self.ts_coords[-1].append(cur_atm_coords)
                            atm_cntr += 1

                line = opened_pwout.readline()
                split_line = None

        # save only the last frame
        if save_all_scf_steps is False:
            self.ts_coords = [self.ts_coords[-1]]

            for key in self.pw_other_info:

                try:
                    self.pw_other_info[key] = [self.pw_other_info[key][-1]]
                except IndexError:
                    pass

    def _convert_frozen_state(self):
        """
        In gaussian frozen atoms are marked using -1, unfrozen atoms are marked using 0.
        In quantum espresso on the other hand, frozen atoms are marked using 0, unfrozen atoms have a 1.
        Since gaussian does not distinguish in which direction atoms are allowed to move, quantum espresso
        does differ between moving directions. Therefor it is necessary to know where the frozen info came
        from and to what it will be converted to.
        #TODO THIS FUNCTION IS JUST A VERY DIRTY FIX UNTIL THE CONCEPT OF
        #TODO HANDLING FROZEN ATOMS IS REDONE (E.G. KNOWING WHERE THE FREEZING INFORMATION COMES FROM)
        """
        # first check if there is a reason to freeze at all
        for atom in self.atoms:
            if atom.ifrz is not None:
                frozen_atom_found = True
                break
        else:
            frozen_atom_found = False

        if frozen_atom_found is True:
            for atom in self.atoms:
                if atom.ifrz == -1:
                    atom.ifrz_x = 0
                    atom.ifrz_y = 0
                    atom.ifrz_z = 0
                elif atom.ifrz == 0 or atom.ifrz is None:
                    atom.ifrz_x = 1
                    atom.ifrz_y = 1
                    atom.ifrz_z = 1
                else:
                    #ifrz is None or 2 for example
                    pass

    def _write_section(self, opened_file_instance, frame_id, keyword):
        """
        Help writing an input section.

        opened_file_instance    file; file to write section to
        frame_id                int; id of coordinates to write
        keyword                 str; pw_entry keyword
        """
        # write entries
        opened_file_instance.write("&{}\n".format(keyword))

        for setting, value in self.pw_entries[keyword].items():

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
        #TODO NUMBER OF ATOMS WRONG WHEN READ FROM GAUSSIAN
        #TODO WRITE FROZEN ATOMS
        """
        self.pw_entries["CONTROL"]["verbosity"] = verbosity
        self.pw_entries["SYSTEM"]["nat"] = len(self.atoms)
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
            for _, atom_type in self.atm_types.items():
                opened_filename.write("{0} {1:>6} {2:}\n".format(
                    atom_type.sitnam, atom_type.weigh,
                    atom_type.pseudopotential))
            opened_filename.write("\n")

            # Atom Coordinates
            opened_filename.write("ATOMIC_POSITIONS {angstrom}\n")
            # check if a single atom is frozen -> all atoms have to be frozen
            for atom in self.atoms:
                if hasattr(atom, "ifrz_x") and hasattr(atom, "ifrz_y") and hasattr(atom, "ifrz_z"):
                    print_frozen_state = True
                    coordinate_string = "{0:<5s} {1:> 18.9f}{2:> 18.9f}{3:> 18.9f}{4:> 18}{5:> 18}{6:> 18}\n"
                    break
            else:
                print_frozen_state = False
                coordinate_string = "{0:<5s} {c[0]:> 18.9f}{c[1]:> 18.9f}{c[2]:> 18.9f}\n"

            for atom, atom_coordinates in zip(self.atoms,
                                              self.ts_coords[frame_id]):

                # write format for frozen atoms
                if print_frozen_state is False:
                    opened_filename.write(coordinate_string.format(
                        atom.sitnam, c=atom_coordinates))
                else:
                    opened_filename.write(coordinate_string.format(
                        atom.sitnam,
                        atom_coordinates[0], atom_coordinates[1], atom_coordinates[2],
                        atom.ifrz_x, atom.ifrz_y, atom.ifrz_z))

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


def read_pw_out(pw_out, read_crystal_sections=False):
    """Read a pw output file."""
    pw_sys = PwStuff()
    pw_sys.read_pwout(pw_out, read_crystal_sections)
    return pw_sys
