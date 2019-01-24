from __future__ import print_function, division
import re
import numpy as np
import scipy.constants
import pdb
import md_stars as mds
import md_universe as mdu
import md_elements as mde
#import log_universe as logu

__version__ = "2018-10-16"

"""
CURRENTLY THIS MODULE IS UNDER CONSTRUCTION AND NOT FULLY FUNCTIONAL!
Forces should also be read.
"""

hartree_eV = scipy.constants.physical_constants["Hartree energy in eV"][0]


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class GauStuff(mdu.Universe):
    """
    Read and write Gaussian files. Totally limited at the moment. Internal
    indices start with 0 (converted during reading and reconverted during
    writing procedure).
    """
    def __init__(self,
                 rwf=None,
                 oldchk=None,
                 chk=None,
                 nproc=None,
                 mem=None,
                 job_settings=None,
                 gaussian_charges=[],
                 gaussian_multiplicities=[],
                 gaussian_other_info={},
                 energy_unit=None):
        """
        > job_settings      str; line with job settings, e.g. #P Opt=tight MP2/6-311++G**
                            or #P SP MP2/6-311++G**
        > coordinate_style       str; cartesian | z-matrix
        >
        >

        Sources:    http://gaussian.com/basissets/
                    http://gaussian.com/dft/
                    http://gaussian.com/mp/
                    http://gaussian.com/capabilities/
                    http://gaussian.com/geom/
                    http://wild.life.nctu.edu.tw/~jsyu/compchem/g09/g09ur/m_molspec.htm

        Note that this is far from complete
        """
        mdu.Universe.__init__(self)

        if rwf is not None:
            self.rwf  =  rwf

        if oldchk is not None:
            self.oldchk = oldchk

        if chk is not None:
            self.chk = chk

        if nproc is not None:
            self.nproc = nproc

        if mem is not None:
            self.mem = mem

        if job_settings is not None:
            self.job_settings = job_settings
        else:
            self.job_settings = ""

        if gaussian_charges != []:
            self.gaussian_charges = gaussian_charges
        else:
            self.gaussian_charges = []

        if gaussian_multiplicities != []:
            self.gaussian_multiplicities = gaussian_multiplicities
        else:
            self.gaussian_multiplicities = []

        if energy_unit is not None:
            self.energy_unit = energy_unit

        if gaussian_other_info != {}:
            self.gaussian_other_info = gaussian_other_info
        else:
            self.gaussian_other_info = {}

    def read_gau(self, gauin, coordinate_style="cartesian", overwrite=False, debug=False):
        """
        """
        self.coordinate_style = coordinate_style

        print("***Gau-Info: Reading  Gaussian-Input-File!")
        cframe = []
        reading = True

        with open(gauin, "r") as gau_in:
            line = gau_in.readline()

            if debug is True:
                print("Parsing link 0 and route section")

            # line != "" -> just to be safe this will eventually finish when
            # loading a wrong file by accident
            while line != "":

                # // LINK 0 SECTION
                if "%oldchk" in line:
                    if not hasattr(self, "oldchk") or overwrite is True:
                        self.oldchk = line.split("=")[1].strip("\n")

                elif "%chk" in line:
                    if not hasattr(self, "chk") or overwrite is True:
                        self.chk = line.split("=")[1].strip("\n")

                elif "%nproc" in line:
                    if not hasattr(self, "nproc") or overwrite is True:
                        self.nproc = int(line.split("=")[1].strip("\n"))

                elif "%mem" in line:
                    if not hasattr(self, "mem") or overwrite is True:
                        self.mem = line.split("=")[1].strip("\n")

                #// ROUTE SECTION (may be scattered over several files)
                elif line.startswith("#"):
                    #pdb.set_trace()

                    # keep reading until empty line is reached, end loop
                    # when it is
                    while line != "\n":
                        if self.job_settings == "" or overwrite is True:
                            self.job_settings += line.rstrip("\n") + " "
                        line = gau_in.readline()

                        # eof reached
                        if line == "":
                            break

                    else:
                        break

                else:
                    pass

                line = gau_in.readline()

            #// TITLE SECTION
            if debug is True:
                print("Title section")

            # skip all empty lines until title line is reached
            while line == "\n":
                line = gau_in.readline()
            else:
                title_line = line
                line = gau_in.readline()

            #// MOLECULE SPECIFICATION SECTION
            if debug is True:
                print("Molecule specification section")

            while line == "\n":
                line = gau_in.readline()
            else:
                charge_mutliplicity_line = re.findall(r'[\w]+', line)

                if self.gaussian_charges == [] or overwrite is True:
                    self.gaussian_charges = [int(i)
                                             for idx, i in enumerate(charge_mutliplicity_line)
                                             if idx % 2 == 0]

                if self.gaussian_multiplicities == [] or overwrite is True:
                    self.gaussian_multiplicities = [int(i)
                                                    for idx, i in enumerate(charge_mutliplicity_line[1:])
                                                    if idx % 2 == 0]

                line = gau_in.readline()

            #// COORDINATES SECTION
            if debug is True:
                print("Coordinates section")

            g_atoms = []
            g_frame = []  # gaussian coordinates
            # counter to substitute given atoms or complement their attributes
            atom_index = 0

            while line != "\n":
                line = line.split()
                columns_coordinates_section_atomic = len(line)

                # element parameters section
                if "Fragment=" in line[0]:
                    fragment_id = re.search(r"Fragment=\d+", line[0]).group(0)
                    fragment_id = int(re.search(r"\d+", fragment_id).group(0))
                    csitnam = line[0].split("(")[0]
                else:
                    fragment_id = None
                    csitnam = line[0]

                # translate atomic number to its element
                if csitnam.isdigit():
                    csitnam = int(csitnam)
                    csitnam = mde.atomicnumber_element[csitnam]

                if "Iso=" in line[0]:
                    print("Iso not implemented yet")

                if "Spin=" in line[0]:
                    print("Iso not implemented yet")

                # element freeze section
                if len(line) > 4 and (line[-4] == "0" or line[-4] == "-1"):
                    cifrz = int(line[-4])
                else:
                    cifrz = None

                # store element info
                cur_atom = mds.Atom(sitnam=csitnam, atm_id=atom_index ,ifrz=cifrz,
                                    grp_id=fragment_id)

                #g_atoms.append(cur_atom)

                try:
                    if overwrite is True:
                        self.atoms[atom_index] = cur_atom
                    else:
                        if not hasattr(self.atoms[atom_index], "sitnam"):
                            self.atoms[atom_index].sitnam = cur_atom.sitnam

                        # add ifrz to current atom if overwrite option is set or ifrz attribute has not been defined yet
                        if hasattr(cur_atom, "ifrz") and not hasattr(self.atoms[atom_index], "ifrz"):
                            self.atoms[atom_index].ifrz = cur_atom.ifrz

                        if hasattr(cur_atom, "grp_id") and not hasattr(self.atoms[atom_index], "grp_id"):
                            self.atoms[atom_index].grp_id = cur_atom.grp_id

                except IndexError:
                    # overwrite the whole entry if one index does not fit
                    self.atoms.append(cur_atom)

                # add coordinates
                if self.coordinate_style == "cartesian":
                    ccoords = np.array([float(i) for i in line[-3:]])
                    g_frame.append(ccoords)
                else:
                    print("Z-Matrix not implemented yet.")

                line = gau_in.readline()
                atom_index += 1

                if line == "":
                    break

            # overwrite or append frame (coordinates) to existing one(s)
            g_frame = np.array(g_frame)

            if self.ts_coords != [] and overwrite is True:
                self.ts_coords[-1] = g_frame
            else:
                self.ts_coords.append(g_frame)

            # skip all empty lines until title line is reached
            while line == "\n":
                line = gau_in.readline()  # first bond line

                # just to be sure this will end
                if line == "":
                    break

            # BOND SECTION
            if debug is True:
                print("Bonds Section")

            if "CONNECTIVITY" in self.job_settings.upper() and bool(re.search(r'(^\d+ \d+ \d+.\d+)|(^\d+\n$)', line)) is True:
                gau_bonds = []

                while line != "\n":
                    line = line.split()

                    # only process line if it really has bond information
                    if len(line) > 1:
                        # get atm-id 1 (always the same for current line)
                        catm_id1 = int(line[0])
                        catm_id1 -= 1  # decrement atom index by 1
                        # read other further bond partners (if present)
                        catms_id2 = [int(i)-1 for idx, i in enumerate(line[1:]) if idx % 2 == 0]
                        # read bond orders (if present)
                        cbnd_orders = [float(i) for idx, i in enumerate(line[2:]) if idx % 2 == 0]

                        # append bond to gaussian given bonds
                        for catm_id2, cbnd_order in zip(catms_id2, cbnd_orders):
                            cbnd = mds.Bond(atm_id1=catm_id1,
                                            atm_id2=catm_id2,
                                            bnd_order=cbnd_order)
                            gau_bonds.append(cbnd)

                    # if there is only one empty line at file bottom
                    line = gau_in.readline()

                    if line == "":
                        break

                # complement attributes that do not exist or overwrite existing ones if wanted
                if len(gau_bonds) == len(self.bonds) and overwrite is False:
                    for idx, gau_bnd in enumerate(gau_bonds):
                        for universe_bnd in self.bonds:
                            if (universe_bnd.atm_id1 == gau_bnd.atm_id1 and
                                universe_bnd.atm_id2 == gau_bnd.atm_id2):
                                # complement attribute or overwrite existing one
                                # if overwriting is active
                                if not hasattr(universe_bnd, "bnd_order"):
                                    universe_bnd.bnd_order = gau_bnd.bnd_order
                                break
                else:
                    self.bonds = gau_bonds

                # create molecules by bond information
                self.fetch_molecules_by_bonds()

            #line = gau_in.readline()

    def write_gau(self, gauout, frame_id, modredundant=None,
                  write_fragments=False, job_settings=None, title=""):
        """
        job_type    modredundant | SP
        method      user choice (e.g. MP2, B3LYP)
                    http://gaussian.com/dft/
                    http://gaussian.com/mp/

        basis_set   user choice (e.g. 6-311+G*)
                    http://gaussian.com/basissets/
        geom        PrintInputOrient|connectivity

        Sources: http://wild.life.nctu.edu.tw/~jsyu/compchem/g09/g09ur/m_molspec.htm
        """
        print("***Gau-Info: Writing Gaussian-Input-File!")

        if job_settings is not None:
            self.job_settings = job_settings

        with open(gauout, "w") as gau_out:
            if hasattr(self, "nproc"):
                gau_out.write("%nproc={}\n".format(self.nproc))

            if hasattr(self, "mem"):
                gau_out.write("%mem={}\n".format(self.mem))

            if hasattr(self, "oldchk"):
                gau_out.write("%oldchk={}\n".format(self.oldchk))

            if hasattr(self, "chk"):
                gau_out.write("%chk={}\n".format(self.chk))

            if hasattr(self, "rwf"):
                gau_out.write("%rwf={}\n".format(self.rwf))

            gau_out.write(self.job_settings + "\n\n")

            # write title line only if allcheck-keyword is not set
            if "ALLCHECK" not in self.job_settings.upper():
                gau_out.write("T: {}\n\n".format(title))

            # write multiplicity and charge line only if allcheck-keyword is not set
            if "ALLCHECK" not in self.job_settings.upper() or (self.charge is not None and self.multiplicity is not None):
                #gau_out.write("{} {}\n".format(self.charge, self.multiplicity))

                for charge, multiplicity in zip(self.gaussian_charges, self.gaussian_multiplicities):
                    gau_out.write("{} {} ".format(charge, multiplicity))

                gau_out.write("\n")

            if self.ts_coords != [] and ("ALLCHECK" not in self.job_settings.upper() or "CHECK" not in self.job_settings.upper()):
                for idx, catm in enumerate(self.atoms):
                    gau_out.write("{}".format(catm.sitnam))

                    if write_fragments is True:
                        for molecule_idx, molecule in enumerate(self.molecules):
                            if idx + 1 in molecule:
                                gau_out.write("(Fragment={})".format(molecule_idx + 1))
                                break

                    # write frozen state of the atom
                    if hasattr(catm, "ifrz") and catm.ifrz is not None:
                        gau_out.write(" {:<3} ".format(catm.ifrz))

                    gau_out.write(" {:> 11.6f}{:> 11.6f}{:> 11.6f}\n".format(
                        self.ts_coords[frame_id][idx][0],
                        self.ts_coords[frame_id][idx][1],
                        self.ts_coords[frame_id][idx][2]))

            gau_out.write("\n")

            # write bonds
            if "CONNECTIVITY" in self.job_settings.upper() and write_fragments is False:
                for cur_atom in self.atoms:
                    gau_out.write("{} ".format(cur_atom.atm_id))
                    # check whether atom is bonded -> write partner atoms and
                    # according bond order
                    for cur_bond in self.bonds:
                        if cur_atom.atm_id == cur_bond.atm_id1:
                            gau_out.write("{} {} ".format(cur_bond.atm_id2,
                                                          cur_bond.bnd_order))
                    gau_out.write("\n")
                gau_out.write("\n")

            if modredundant is not None:
                gau_out.write(modredundant)
                gau_out.write("\n" * 2)

            gau_out.write("\n" * 4)

    def _read_gau_log_summary(self, result_str, overwrite=False):
        """
        """
        result_str = result_str.split("\\\\")
        atoms_coords = result_str[3].split("\\")
        # remove first entry which is charge and multiplicity
        atoms_coords.pop(0)
        ts_coords = []

        for atom_coords in atoms_coords:
            atom_coords = atom_coords.split(",")
            coords = np.array([float(i) for i in atom_coords[-3:]])
            ts_coords.append(coords)

        # get number of imaginary frequencies (if frequencies were calculated)
        for subresult in result_str:
            if "NImag" in subresult:
                self.gaussian_other_info["NImag"] = int(subresult.split("NImag=")[1])

        if overwrite is False:
            self.ts_coords.append(ts_coords)
        else:
            self.ts_coords[-1] = ts_coords

    def read_gau_log(self, gau_log, save_all_scf_steps=False, overwrite=False, read_summary=False):
        """
        Read the last coordinates from a gaussian log file.

        Overwrite overwrites the last frame.
        """
        #TODO read atom by initial coordinates and not by scf cycles

        if overwrite is True:
            self.ts_coords = []

        #print("Reading last frame of the output file.")
        if "scf_energies" not in self.gaussian_other_info or overwrite is True:
            self.gaussian_other_info["scf_energies"] = []

        all_scf_cycles_coords = []
        current_scf_cycles_coords = []
        all_scf_cycles_energies = []
        scf_cycles_energies = []
        log_resume = ""
        g_atoms = []
        atom_index = 0
        read_element_numbers = True
        #scanned_coordinates = []

        # read geometries from all scf cycles and their corresponding energies
        with open(gau_log, "r") as opened_gau_log:
            line = opened_gau_log.readline()

            # read all scf cycles and the corresponding energy
            while line != "":

                if "Standard orientation" in line:
                    cframe = []
                    # skip the following 4 lines
                    for _ in xrange(5):
                        line = opened_gau_log.readline()

                    while not line.startswith(" ----------------------------"):
                        split_line = line.split()
                        coords = [float(i) for i in split_line[3:]]
                        coords = np.array(coords)
                        cframe.append(coords)

                        # get element number, i.e. element name
                        if read_element_numbers is True:
                            element_number = int(split_line[1])
                            element_name = mde.element_name[element_number]
                            cur_atom = mds.Atom(sitnam=element_name,
                                                atm_id=atom_index)
                            g_atoms.append(cur_atom)
                            atom_index += 1

                        line = opened_gau_log.readline()

                    # append current frame to current scf_cycle
                    current_scf_cycles_coords.append(cframe)
                    # stop reading the element numbers after first entry
                    read_element_numbers = False

                elif "SCF Done:" in line:
                    scf_energy = hartree_eV * float(line.split()[4])
                    scf_cycles_energies.append(scf_energy)

                elif "!   Optimized Parameters   !" in line:
                    # current entry is finished
                    all_scf_cycles_coords.append(current_scf_cycles_coords)
                    all_scf_cycles_energies.append(scf_cycles_energies)
                    # reset current optimized scf cycles
                    current_scf_cycles_coords = []
                    # reset energies for next run to read
                    scf_cycles_energies = []
                    #self.gaussian_other_info("optimized_parameters")

                    # skip the next 5 lines
                    for _ in xrange(6):
                        line = opened_gau_log.readline()

                    while not line.startswith(" ----------------------------"):
                        split_line = line.split()

                        parameter_definition = split_line[2]
                        parameter_value = float(split_line[3])
                        self.gaussian_other_info[parameter_definition].append(parameter_value)

                        #for scanned_coordinate in scanned_coordinates:
                        #    if split_line[2] == scanned_coordinate:
                        #        self.gaussian_other_info[scanned_coordinate].append(float(split_line[3]))

                        line = opened_gau_log.readline()

                # get scanned coordinates (if mod redundant is used)
                elif "Initial Parameters" in line:
                    # skip the next 4 lines
                    for _ in xrange(5):
                        line = opened_gau_log.readline()

                    while not line.startswith(" ----------------------------"):
                        split_line = line.split()

                        # prepare containers for current parameter definition
                        parameter_definition = split_line[2]
                        self.gaussian_other_info[parameter_definition] = []

                        #if "Scan" in line:
                        #    scanned_coordinates.append(split_line[2])

                        line = opened_gau_log.readline()

                    #for scanned_coordinate in scanned_coordinates:
                    #    self.gaussian_other_info[scanned_coordinate] = []

                # read the summary
                elif line.startswith(" 1\\1\\") is True:

                    while not line.endswith("@\n"):
                        line = line.lstrip()
                        line = line.rstrip()
                        # create one huge string since lines could be oddly wrapped
                        log_resume += line
                        line = opened_gau_log.readline()
                else:
                    pass

                line = opened_gau_log.readline()

        #return (all_scf_cycles_coords, all_scf_cycles_energies)
        # extract geometries with lowest energy, omit other scf geometries
        # and energies
        for scf_cycle_coords, scf_cycle_energy in zip(all_scf_cycles_coords,
                                                      all_scf_cycles_energies):
            cur_cycle_min_scf_energy = 1e12
            cur_cycle_min_scf_energy_idx = None

            # get index of frame with lowest energy
            for index_energy, cur_energy in enumerate(scf_cycle_energy):
                if cur_energy < cur_cycle_min_scf_energy:
                    cur_cycle_min_scf_energy = cur_energy
                    cur_cycle_min_scf_energy_idx = index_energy

            try:
                self.ts_coords.append(scf_cycle_coords[cur_cycle_min_scf_energy_idx])
                self.gaussian_other_info["scf_energies"].append(scf_cycle_energy[cur_cycle_min_scf_energy_idx])
            except IndexError:
                print("***Warning: Gaussian run was aborted before it finished")

        if read_summary is True:
            self._read_gau_log_summary(log_resume, overwrite=overwrite)

        # # split string by its entries
        # log_resume = log_resume.split("\\")
        # del line
        # del get_lines

        # g_atoms = []
        # g_frame = []
        # # switch which defines that coordinates section is reached
        # coordinates_cntr = 0
        # # atom indx
        # cidx = 0

        # for entry in log_resume:
        #    if entry == "":
        #        coordinates_cntr += 1

        #    if coordinates_cntr == 3:
        #        entry = entry.split(",")

        #        if len(entry) == 4 or len(entry) == 5:
        #            if len(entry) == 4:
        #                ccoords = np.array([float(i) for i in entry[1:]])
        #            else:
        #                ccoords = np.array([float(i) for i in entry[2:]])

        #            csitnam = entry[0]

        #            # store info
        #            cur_atom = mds.Atom(sitnam=csitnam, atm_id=cidx)

        #            g_atoms.append(cur_atom)
        #            g_frame.append(ccoords)

        #            # increment atom index
        #            cidx += 1

        #    if "NImag" in entry:
        #        # get number of imaginary frequencies
        #        self.gaussian_other_info["NImag"] = int(entry.split("NImag=")[1])

        #    if "Version" in entry:
        #        other_entries = entry.split("\\")
        #        for other_entry in other_entries:
        #            if "HF" in other_entry:
        #                energies_entry = other_entry.split(",")
        #                energies_entry = [float(i)*hartree_eV for i in energies_entry]
        #                self.gaussian_other_info["energies_entry"] = energies_entry
#
        #del cidx
        #del coordinates_cntr
#
        #g_frame = np.array(g_frame)
#
        for idx, gatm in enumerate(g_atoms):
            try:
                if overwrite is True:
                    self.atoms[idx] = gatm
                else:
                    if not hasattr(self.atoms[idx], "sitnam"):
                        self.atoms[idx].sitnam = gatm.sitnam

                    if not hasattr(self.atoms[idx], "atm_id"):
                        self.atoms[idx].atm_id = gatm.atm_id

            except IndexError:
                # overwrite the whole entry if one index does not fit
                self.atoms = g_atoms
                break

        #if self.ts_coords != [] and overwrite is True:
        #    self.ts_coords[-1] = g_frame
        #else:
        #    self.ts_coords.append(g_frame)


def read_gauin(gauin, coordinate_style="cartesian", overwrite=False, debug=False):
    """
    Read the gaussian input file.

    Parameters
    ----------
    gauin : str
        name of the gaussian input file

    Returns
    -------
    gau_sys : GauStuff object
        An object of GauStuff which can be further processed

    """
    gau_sys = GauStuff()
    gau_sys.read_gau(gauin, coordinate_style, overwrite, debug)
    return gau_sys
