from __future__ import print_function, division
import os
import re
import numpy as np
import scipy.constants as sc
import time
import md_box as mdb
import md_stars as mds
import md_universe as mdu

__version__ = "2018-05-03"

bohr_angstrom = sc.value("Bohr radius")/sc.angstrom
angstrom_bohr = sc.angstrom/sc.value("Bohr radius")

class PwStuff(mdu.Universe):
    """
    PW Stuff
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
                        if re.match("^\d+?\.\d+?$", split_line[1]):
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
                        if re.match("^\d+?\.\d+?$", split_line[1]):
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
                        if re.match("^\d+?\.\d+?$", split_line[1]):
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
                        if re.match("^\d+?\.\d+?$", split_line[1]):
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
                        if re.match("^\d+?\.\d+?$", split_line[1]):
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
                    self.pw_entries["K_POINTS"]["k_point_grid"] = []

                    while line != '':
                        line = opened_pwin.readline()

                        if line.startswith("\n"):
                            break

                        split_line = line.split()

                        if kpoints_line[1] == "automatic":
                            self.pw_entries["K_POINTS"]["k_point_grid"].append([int(i) for i in split_line])
                        elif kpoints_line[1] == "gamma":
                            self.pw_entries["K_POINTS"]["k_point_grid"] = ["\n"]
                            break
                        else:
                            pass

                elif line.startswith("CELL_PARAMETERS"):
                    # get box vectors
                    cbox = mdb.Box(boxtype="cartesian")

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
                    del line_cntr
                else:
                    pass

                line = opened_pwin.readline()

        # convert celldm (= alat) to box vector a with angstrom
        try:
            self.ts_boxes[-1].ltc_a = float(self.pw_entries["SYSTEM"]["celldm(1)"]*bohr_angstrom)
        except KeyError:
            pass

        # final check
        if len(self.atoms) != self.pw_entries["SYSTEM"]["nat"]:
            print("***Warning: Number of atoms and SYSTEM entry 'nat' differ!")
            #time.sleep(5)

        if os.path.isdir("/home/hpc/bccc/bccc34/opt/Programs/quantum_espresso_pseudo_potentials") is False:
            print("***Warning: Folder for Pseudopotentials does not exist!")
            #time.sleep(5)

    def read_pwout(self, pwout):
        """
        CAVEAT: UNDER CONSTRUCTION! Read the output of pw.x
        Currently this only reads the coordinates and cell vectors.
        Cell vector alat (celldm(1)) is converted to angstrom when read.
        """
        with open(pwout) as opened_pwout:
            line = opened_pwout.readline()
            while line != '':
                #TODO box vectors should be converted to lattice or their real lengths
                if line.startswith("CELL_PARAMETERS"):
                    # get alat
                    split_line = line.split()

                    #TODO calculate real cartesian box vectors from given ones
                    # get box vectors
                    cbox = mdb.Box(boxtype="cartesian")

                    if "alat" in line:
                        cbox.ltc_a = float(split_line[2].strip(")"))*bohr_angstrom
                        cbox.crt_a = [float(i) for i in opened_pwout.readline().split()]
                        cbox.crt_b = [float(i) for i in opened_pwout.readline().split()]
                        cbox.crt_c = [float(i) for i in opened_pwout.readline().split()]

                    self.ts_boxes.append(cbox)
                elif line.startswith("ATOMIC_POSITIONS"):
                    # prepare container for coordinates to come
                    self.ts_coords.append([])

                    # read the coordinates
                    while line != '':
                        line = opened_pwout.readline()

                        if line.startswith("\n") or line.startswith("End final coordinates"):
                            break

                        split_line = line.split()

                        cur_atm = mds.Atom(sitnam=split_line[0])
                        cur_atm_coords = np.array([float(i) for i in split_line[1:]])
                        self.ts_coords[-1].append(cur_atm_coords)
                else:
                    pass

                line = opened_pwout.readline()

    def write_pwin(self, frame_id, filename):
        """
        Write an input file for pw.x. alat is converted to angstrom when read
        """
        with open(filename, "w") as opened_filename:

            # write entries
            for keyword in self.pw_entries:

                # skip K_POINTS entry since it will be used later on
                if keyword is "K_POINTS":
                    continue
                else:
                    opened_filename.write("&{}\n".format(keyword))

                    for setting, value in self.pw_entries[keyword].iteritems():

                        # calculate celldm(1), which derives from lattice vector a
                        if setting == "celldm(1)":
                            opened_filename.write("    {0:<24s}= {1}\n".format(setting, self.ts_boxes[frame_id].ltc_a*angstrom_bohr))
                        else:
                            opened_filename.write("    {0:<24s}= {1}\n".format(setting, value))

                    opened_filename.write("/\n\n")

            # Pseudopotentials (atom types)
            opened_filename.write("ATOMIC_SPECIES\n")
            for index, atom_type in self.atm_types.iteritems():
                opened_filename.write("{0} {1:>6} {2:}\n".format(atom_type.sitnam,
                                                               atom_type.weigh,
                                                               atom_type.pseudopotential))
            opened_filename.write("\n")

            # Atom Coordinates
            #TODO currently atomic positions are only supported in angstrom
            opened_filename.write("ATOMIC_POSITIONS {angstrom}\n")
            for atom, atom_coordinates in zip(self.atoms, self.ts_coords[frame_id]):
                opened_filename.write("{0:<5s} {c[0]:> 18.9f}{c[1]:> 18.9f}{c[2]:> 18.9f}\n".format(atom.sitnam, c=atom_coordinates))
            opened_filename.write("\n")

            # K Points
            opened_filename.write("K_POINTS {}\n".format(self.pw_entries["K_POINTS"]["option"]))

            # write grid or skip entry if gamma point is set
            if self.pw_entries["K_POINTS"]["option"] == "automatic":
                opened_filename.write("{}\n".format(self.pw_entries["K_POINTS"]["k_point_grid"]))

            opened_filename.write("\n")

            # Cell
            #TODO currently only alat as cell option is supported
            opened_filename.write("CELL_PARAMETERS {alat}\n")
            opened_filename.write("{v[0]:> 24.9f}{v[1]:> 18.9f}{v[2]:> 18.9f}\n".format(v=self.ts_boxes[frame_id].crt_a))
            opened_filename.write("{v[0]:> 24.9f}{v[1]:> 18.9f}{v[2]:> 18.9f}\n".format(v=self.ts_boxes[frame_id].crt_b))
            opened_filename.write("{v[0]:> 24.9f}{v[1]:> 18.9f}{v[2]:> 18.9f}\n".format(v=self.ts_boxes[frame_id].crt_c))
            opened_filename.write("\n")
