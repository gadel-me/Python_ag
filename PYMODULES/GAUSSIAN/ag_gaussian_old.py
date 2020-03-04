import re
import numpy as np
import md_stars as mds
import md_universe as mdu


class GauStuff(mdu.Universe):
    """
    Read and write Gaussian files. Totally limited at the moment. Internal
    indices start with 0 (converted during reading and reconverted during
    writing procedure).
    """

    def __init__(
        self,
        chk=None,
        nproc=None,
        mem=None,
        job_type=None,
        method=None,
        basis_set=None,
        geom=None,
        charge=None,
        dispersion=None,
        multiplicity=None,
        energy_unit=None,
    ):
        """
        Sources:    http://gaussian.com/basissets/
                    http://gaussian.com/dft/
                    http://gaussian.com/mp/
                    http://gaussian.com/capabilities/
                    http://gaussian.com/geom/

        Note that this is far from complete
        """
        mdu.Universe.__init__(self)

        if chk is not None:
            self.chk = chk
        if nproc is not None:
            self.nproc = nproc
        if mem is not None:
            self.mem = mem
        if job_type is not None:
            self.job_type = job_type
        if method is not None:
            self.method = method
        if basis_set is not None:
            self.basis_set = basis_set
        if geom is not None:
            self.geom = geom
        if charge is not None:
            self.charge = charge
        if multiplicity is not None:
            self.multiplicity = multiplicity
        if energy_unit is not None:
            self.energy_unit = energy_unit
        if dispersion is not None:
            self.dispersion = dispersion

        self.ab_initio_keywords = set(
            (
                "HF",
                "RHF",
                "UHF",
                "ROHF",
                "OSS",
                "GVB",
                "CASSCF",
                "MP2",
                "MP3",
                "MP4",
                "MP4DQ",
                "MP$SDQ",
                "MP4SDTQ",
                "CI",
                "CIS",
                "CISD",
                "QCISD",
                "QCISD(T)",
            )
        )
        self.basis_sets = set(
            (
                "STO-3G",
                "3-21G",
                "4-21G",
                "6-21G",
                "6-31G",
                "LP-31G",
                "LP-41G",
                "6-311G",
                "MC-311G",
                "D95",
                "D95V",
                "SEC",
                "CEP-4G",
                "CEP-31G",
                "CEP-121G",
                "LANLIMB",
                "LANLIDZ",
            )
        )
        self.job_types = set(("SP", "Opt", "Freq"))
        self.gau_energies = []

    def read_gau(self, gauin, overwrite=False):
        """
        """
        print("***Gau-Info: Reading  Gaussian-Input-File!")

        (chk, nproc, mem, job_type, method, basis_set, geom, charge, multiplicity) = [
            None for _ in range(9)
        ]
        cframe = []

        with open(gauin, "r") as gau_in:
            for line in gau_in:
                if "%chk" in line:
                    if not hasattr(self, "chk"):
                        self.chk = line.split("=")[1].strip("\n")
                elif "%nproc" in line:
                    if not hasattr(self, "nproc"):
                        self.nproc = int(line.split("=")[1].strip("\n"))
                elif "%mem" in line:
                    if not hasattr(self, "mem"):
                        mem = line.split("=")[1].strip("\n")
                        self.mem = int(re.findall(r"^\d+", mem)[0])
                elif line.startswith("#"):
                    # job settings direction
                    job_settings = line.split()

                    # classify job settings
                    for setting in job_settings:
                        # remove # e.g. #SP
                        if setting.startswith("#"):
                            setting = setting.strip("#")

                        if setting in self.job_types:
                            if not hasattr(self, "job_type"):
                                self.job_type = setting
                        elif "geom" in setting:
                            if not hasattr(self, "geom"):
                                self.geom = setting.split("=")[1]
                        elif "EmpiricalDispersion" in setting:
                            if not hasattr(self, "dispersion"):
                                self.dispersion = setting.split("=")[1]
                        elif "/" in setting:
                            setting = setting.split("/")
                            if not hasattr(self, "method"):
                                self.method = setting[0]
                            if not hasattr(self, "basis_set"):
                                self.basis_set = setting[1]
                        else:
                            pass

                elif re.findall(r"^-?\d+ \d$", line):
                    line = line.split()

                    if not hasattr(self, "charge"):
                        self.charge = int(line[0])

                    if not hasattr(self, "multiplicity"):
                        self.multiplicity = int(line[1])

                    read_coords_section = True

                    while read_coords_section is True:

                        # read info
                        try:
                            line = next(gau_in)
                        except StopIteration:
                            # last line of document reached
                            break

                        # quit reading section when empty line occurs
                        if line == "\n":
                            read_coords_section = False
                            break

                        line = line.split()
                        csitnam = line[0]
                        ccoords = np.array([float(i) for i in line[1:]])

                        # store info
                        cur_atom = mds.Atom(sitnam=csitnam)
                        self.atoms.append(cur_atom)
                        cframe.append(ccoords)

                    self.ts_coords.append(cframe)

                elif hasattr(self, "geom") and self.geom == "connectivity":
                    # section with bonds-information
                    if bool(re.search(r"(^\d+ \d+ \d+.\d+)|(^\d+\n$)", line)) is True:
                        line = line.split()

                        # get rest of the bond-entries
                        read_connectivity_section = True

                        while read_connectivity_section is True:

                            # stop reading current entry if line is empty
                            if line == []:
                                break

                            # only process line if it really has bond information
                            if len(line) > 1:
                                # get atm-id 1 (always the same for current line)
                                catm_id1 = int(line[0])

                                for idx, cur_subentry in enumerate(line[1:]):

                                    # get atm-id 2
                                    if idx % 2 == 0:
                                        catm_id2 = int(cur_subentry)
                                        # decrease indices by 1 since internally
                                        # atom-indices start with 0
                                        cbnd = mds.Bond(
                                            atm_id1=catm_id1 - 1, atm_id2=catm_id2 - 1
                                        )
                                    else:
                                        cbnd_order = float(cur_subentry)
                                        cbnd.bnd_order = cbnd_order
                                        self.bonds.append(cbnd)

                            # if there is only one empty line at file bottom
                            try:
                                line = gau_in.next().split()
                            except StopIteration:
                                break

                else:
                    pass

    def write_gau(self, gauout, frame_id, title=""):
        """
        job_type    modredundant | SP
        method      user choice (e.g. MP2, B3LYP)
                    http://gaussian.com/dft/
                    http://gaussian.com/mp/

        basis_set   user choice (e.g. 6-311+G*)
                    http://gaussian.com/basissets/
        geom        PrintInputOrient|connectivity
        """
        print("***Gau-Info: Writing Gaussian-Input-File!")

        with open(gauout, "w") as gau_out:
            gau_out.write("%chk={}\n".format(self.chk))
            gau_out.write("%nproc={}\n".format(self.nproc))
            gau_out.write("%mem={}GB\n".format(self.mem))

            if hasattr(self, "geom"):
                gau_out.write(
                    "#P {} {}/{} geom={} ".format(
                        self.job_type, self.method, self.basis_set, self.geom
                    )
                )
            else:
                gau_out.write(
                    "#P {} {}/{} ".format(self.job_type, self.method, self.basis_set)
                )

            if hasattr(self, "dispersion"):
                gau_out.write("EmpiricalDispersion={}".format(self.dispersion))

            gau_out.write("\n")

            gau_out.write("\n")
            gau_out.write("Title: {}\n".format(title))
            gau_out.write("\n")
            gau_out.write("{} {}\n".format(self.charge, self.multiplicity))

            for idx, catm in enumerate(self.atoms):
                gau_out.write(
                    "{:<3} {:> 11.6f}{:> 11.6f}{:> 11.6f}\n".format(
                        catm.sitnam,
                        self.ts_coords[frame_id][idx][0],
                        self.ts_coords[frame_id][idx][1],
                        self.ts_coords[frame_id][idx][2],
                    )
                )

            gau_out.write("\n")

            # write bonds
            if hasattr(self, "geom") and self.geom.upper() == "CONNECTIVITY":
                # since gaussian wants bonds in a kind of condensed style,
                # prepare a dictionary to do so
                a = list(range(len(self.atoms)))
                b = [[] for _ in range(len(self.atoms))]
                d = dict(list(zip(a, b)))

                # merge bonds of same atom 1 utilizing a dictionary
                for i in range(len(self.bonds)):
                    d[self.bonds[i].atm_id1].append(self.bonds[i])

                # write bond-entry
                for i in d:
                    # write current atom
                    gau_out.write("{}".format(i + 1))

                    for val in d[i]:
                        if val != []:
                            gau_out.write(
                                " {} {}".format(val.atm_id2 + 1, float(val.bnd_order))
                            )

                    # newline after each entry
                    gau_out.write("\n")

    def read_gau_log(self, gau_log):
        """
        Read the last coordinates from a gaussian log file.
        """
        log_resume = ""
        get_lines = False
        with open(gau_log, "r") as gau_log:
            for line in gau_log:
                if line.startswith(" 1\\1\\") is True:
                    get_lines = True
                elif line.endswith("\\@\n") is True:
                    get_lines = False
                else:
                    pass

                if get_lines is True:
                    line = line.lstrip()
                    line = line.rstrip()
                    # create one huge string since lines could be oddly wrapped
                    log_resume += line

        # split string by its entries
        log_resume = log_resume.split("\\")
        del line
        del get_lines

        cframe = []
        coordinates_cntr = 0

        for entry in log_resume:
            if entry == "":
                coordinates_cntr += 1

            if coordinates_cntr == 3:
                entry = entry.split(",")
                if len(entry) == 4:
                    csitnam = entry[0]
                    ccoords = np.array([float(i) for i in entry[1:4]])

                    # store info
                    cur_atom = mds.Atom(sitnam=csitnam)
                    self.atoms.append(cur_atom)
                    cframe.append(ccoords)

        self.ts_coords.append(cframe)
