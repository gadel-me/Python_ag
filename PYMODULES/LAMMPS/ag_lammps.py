from __future__ import print_function
import pdb
import os
import numpy as np
import copy
import md_box as mdb
import md_stars as mds
import md_universe as mdu
import ag_vectalg as agv
import struct
import ag_lmpdcd_helpers as agldh
from natsort import natsorted
#import collections

__version__ = "2018-10-25"


class LmpStuff(mdu.Universe):
    """
    Lammps-Data-File-Class
    """
    def __init__(self):
        """
        Initialize general stuff.
        """
        # generates lists for atom-, bond-, angle-types, etc. pp.
        mdu.Universe.__init__(self)

    def read_lmpdat(self, lmpdat, energy_unit=None, angle_unit=None,
                    overwrite_data=False, debug=False):
        """
        energy_unit eV, kCal, kJ
        angle_unit  deg, rad
        cgcmm       boolean; parse cgcmm although not given in file header;
                    bond coefficients may be given by numbers 1, 2, 3 in
                    comment line (first entry) after each entry in the Bond Coeffs section:
                    1 == single bond
                    2 == double bond
                    3 == triple bond
                    the atom types that are connected by that bond may also be given
                    by the 2nd and 3rd element in the comment (but are currently not read)
                    atom types may be declared by a comment in the Masses
                    section
                    atom names and molecule names may be given by a comment in
                    the Atoms section, e.g.:
                    1 1 1 -0.141306 -3.686300 -1.052600 0.749600 # C1 cbz
        """
        pair_ii = False

        # check if there are frames existing before data is loaded
        if self.ts_coords:
            print("***Info Loading coordinates from data-file on top of " +
                  "already loaded coordinates!")

        lmpdat_box = mdb.Box(boxtype="lammps")

        with open(lmpdat, "r") as lmpdat_in:
            line = lmpdat_in.readline()

            if "CGCMM" in line.upper() and debug is True:
                print("***Info: CGCMM-Style found! " +
                      "Trying to parse additional data.")

            for line in lmpdat_in:
                # /// general stuff ///
                if "atoms" in line:
                    total_atms = int(line.split()[0])
                elif "bonds" in line:
                    total_bnds = int(line.split()[0])
                elif "angles" in line:
                    total_angs = int(line.split()[0])
                elif "dihedrals" in line:
                    total_dihs = int(line.split()[0])
                elif "impropers" in line:
                    total_imps = int(line.split()[0])
                elif "atom types" in line:
                    total_atmtypes = int(line.split()[0])
                elif "bond types" in line:
                    total_bndtypes = int(line.split()[0])
                elif "angle types" in line:
                    total_angtypes = int(line.split()[0])
                elif "dihedral types" in line:
                    total_dihtypes = int(line.split()[0])
                elif "improper types" in line:
                    total_imptypes = int(line.split()[0])

                # /// box settings ///
                elif "xlo xhi" in line:
                    line = line.split()
                    lmpdat_box.lmp_xlo = float(line[0])
                    lmpdat_box.lmp_xhi = float(line[1])
                elif "ylo yhi" in line:
                    line = line.split()
                    lmpdat_box.lmp_ylo = float(line[0])
                    lmpdat_box.lmp_yhi = float(line[1])
                elif "zlo zhi" in line:
                    line = line.split()
                    lmpdat_box.lmp_zlo = float(line[0])
                    lmpdat_box.lmp_zhi = float(line[1])
                elif "xy xz yz" in line:
                    line = line.split()
                    lmpdat_box.lmp_xy = float(line[0])
                    lmpdat_box.lmp_xz = float(line[1])
                    lmpdat_box.lmp_yz = float(line[2])

                # /// atom types (masses) ///
                elif "Masses" in line:
                    lmpdat_in.next()  # skip empty line

                    # parse mass entry
                    atm_tp_old_new = {}
                    for atmcnt in xrange(total_atmtypes):
                        line = lmpdat_in.next()
                        cur_atype = mds.Atom()
                        # parse cgcmm-section
                        lmpdat_stuff, csitnam, cres = self._parse_cgcmm(line)

                        if csitnam is not None:
                            cur_atype.sitnam = csitnam

                        atm_key = int(lmpdat_stuff[0])
                        cur_atype.weigh = float(lmpdat_stuff[1])
                        self.atm_types[atmcnt] = cur_atype
                        atm_tp_old_new[atm_key] = atmcnt

                # /// bond types(coeffs) ///
                elif "Bond Coeffs" in line:
                    lmpdat_in.next()  # skip empty line

                    # parse bond-type entries
                    bnd_tp_old_new = {}
                    for bndcnt in xrange(total_bndtypes):
                        line = lmpdat_in.next()
                        line, comment = self._split_line(line)  # split line into data and comment
                        line = line.split()

                        cur_bndtype = mds.Bond()

                        if comment is not None:
                            cur_bndtype.comment = comment.rstrip()

                        cur_bndtype.energy_unit = energy_unit
                        bnd_key = int(line[0])

                        if len(line) > 2:
                            cur_bndtype.prm1 = float(line[1])
                            cur_bndtype.prm2 = float(line[2])

                        if len(line) > 3:
                            cur_bndtype.prm3 = float(line[3])

                        if len(line) > 4:
                            cur_bndtype.prm4 = float(line[4])

                        # read bond order if given
                        if comment is not None:
                            comment = comment.split()
                            # read bond order if given, skip other comments
                            try:
                                cur_bndtype.bnd_order = int(comment[0])
                            except ValueError:
                                pass

                            # read atom types the bond is between
                            if len(comment) > 2:
                                cur_bndtype.sitnam_2 = comment[2]
                            if len(comment) > 1:
                                cur_bndtype.sitnam_1 = comment[1]

                        self.bnd_types[bndcnt] = cur_bndtype
                        bnd_tp_old_new[bnd_key] = bndcnt

                        # Check if force field units are o.k.
                        if debug is True:
                            self.bnd_types[bndcnt].check_bnd_type()

                # /// angle types(coeffs) ///
                elif "Angle Coeffs" in line:
                    lmpdat_in.next()  # skip empty line

                    # parse angle-type entries
                    ang_tp_old_new = {}
                    for angcnt in xrange(total_angtypes):
                        line = lmpdat_in.next()
                        line, comment = self._split_line(line)
                        line = line.split()

                        cur_angtype = mds.Angle()

                        if comment is not None:
                            cur_angtype.comment = comment.rstrip()

                        cur_angtype.energy_unit = energy_unit
                        cur_angtype.angle_unit = angle_unit
                        ang_key = int(line[0])
                        #print(ang_key)

                        if len(line) > 2:
                            cur_angtype.prm1 = float(line[1])
                            cur_angtype.prm2 = float(line[2])

                        if len(line) > 3:
                            cur_angtype.prm3 = float(line[3])

                        if len(line) > 4:
                            cur_angtype.prm4 = float(line[4])

                        self.ang_types[angcnt] = cur_angtype
                        ang_tp_old_new[ang_key] = angcnt

                        # check force field units
                        self.ang_types[angcnt].check_ang_type()

                # /// dihedral types(coeffs) ///
                elif "Dihedral Coeffs" in line:
                    if debug is True:
                        print("***Lammps-Data-Info: Only charmm-dihedral-style supported (atm)!")
                    lmpdat_in.next()  # skip empty line

                    # parse dihedral-type entries
                    dih_tp_old_new = {}
                    for dihcnt in xrange(total_dihtypes):
                        line = lmpdat_in.next()
                        line, comment = self._split_line(line)
                        line = line.split()

                        cur_dihtype = mds.Dihedral()

                        if comment is not None:
                            cur_dihtype.comment = comment.rstrip()

                        cur_dihtype.energy_unit = energy_unit
                        dih_key = int(line[0])
                        cur_dihtype.prm_k = float(line[1])
                        cur_dihtype.prm_n = int(line[2])
                        cur_dihtype.prm_d = float(line[3])
                        cur_dihtype.weigh_factor = float(line[4])

                        self.dih_types[dihcnt] = cur_dihtype
                        dih_tp_old_new[dih_key] = dihcnt

                        # check force field units
                        if debug is True:
                            self.dih_types[dihcnt].check_dih_type()

                # /// improper types(coeffs) ///
                elif "Improper Coeffs" in line:
                    if debug is True:
                        print("***Lammps-Data-Info: Only cvff-improper-style supported (atm)!")

                    lmpdat_in.next()

                    # parse improper-type entries
                    imp_tp_old_new = {}
                    for impcnt in xrange(total_imptypes):
                        line = lmpdat_in.next()
                        line, comment = self._split_line(line)
                        line = line.split()

                        cur_imptype = mds.Improper()

                        if comment is not None:
                            cur_imptype.comment = comment.rstrip()

                        cur_imptype.energy_unit = energy_unit
                        imp_key = int(line[0])
                        cur_imptype.prm_k = float(line[1])
                        cur_imptype.prm_d = int(line[2])
                        cur_imptype.prm_n = int(line[3])

                        self.imp_types[impcnt] = cur_imptype
                        imp_tp_old_new[imp_key] = impcnt

                        # check force field units
                        if debug is True:
                            self.imp_types[impcnt].check_imp_type()

                # /// pair coefficients entry ///
                elif "Pair Coeffs" in line:
                    if debug is True:
                        print("***Lammps-Data-Info: Parsing Pair Coeffs")

                    lmpdat_in.next()
                    total_pairtypes = total_atmtypes

                    for _ in xrange(total_pairtypes):
                        line = lmpdat_in.next()
                        line, comment = self._split_line(line)
                        line = line.split()

                        # atm type, epsilon ii, sigma ii
                        atm_key_i  = atm_tp_old_new[int(line[0])]
                        epsilon_ii = float(line[1])
                        sigma_ii   = float(line[2])

                        #cur_pairtype = mds.LongRange(
                        #    atm_key_i=atm_key_i,  # old atom-id to new one
                        #    epsilon_ij=epsilon_ii,
                        #    sigma_ij=sigma_ii,
                        #    pairs="ii"
                        #)
                        #self.pair_types.append(cur_pairtype)

                        # also give info to atm_types; decrease atom-index by 1
                        # since internally we are starting at 0
                        self.atm_types[atm_key_i].epsilon = epsilon_ii
                        self.atm_types[atm_key_i].sigma   = sigma_ii

                    pair_ii = True

                elif "PairIJ Coeffs" in line:
                    if debug is True:
                        print("***Lmpdat-Info: Parsing PairIJ Coeffs")

                    lmpdat_in.next()
                    total_pairtypes = total_atmtypes*(total_atmtypes+1)/2

                    for _ in xrange(total_pairtypes):
                        line = lmpdat_in.next()
                        line, comment = self._split_line(line)
                        line = line.split()
                        atm_key_i  = atm_tp_old_new[int(line[0])]
                        atm_key_j  = atm_tp_old_new[int(line[1])]
                        epsilon_ij = float(line[2])
                        sigma_ij   = float(line[3])

                        # Still read IJ pair coeffs since they might differ from
                        # IJ mixing procedures
                        cur_pairtype = mds.LongRange(
                            atm_key_i=atm_key_i,  # convert atom id from input to internal atom id
                            atm_key_j=atm_key_j,
                            epsilon_ij=epsilon_ij,
                            sigma_ij=sigma_ij,
                            pairs="ij"
                        )
                        self.pair_types.append(cur_pairtype)

                        # get pair coefficients for each atom type (i.e. i == j)
                        if atm_key_i == atm_key_j:
                            self.atm_types[atm_key_i].epsilon = epsilon_ij
                            self.atm_types[atm_key_i].sigma   = sigma_ij

                # /// atoms entry ///
                elif "Atoms" in line:
                    # read whole section first, sort by id, then re-index the atom-ids
                    #self.atm_idx_id = {}
                    atm_id_old_new = {}
                    tmp_ts = []
                    lmpdat_in.next()  # skip empty line
                    tmp_atm_lines = []

                    # parse 'Atoms' and append to self.atoms
                    for atmcnt in xrange(total_atms):
                        line = lmpdat_in.next()

                        # divide line in lmpdat- and cgcmm-stuff
                        #lmpdat_stuff, csitnam, cres = self._parse_cgcmm(line)  #TODO old garbage
                        tmp_atm_lines.append(self._parse_cgcmm(line))

                    # sort lines by id
                    if debug is True:
                        print("***Lammps-Data-Info: Sorting atoms by their id's, " +
                              "starting with the smallest one from given data file.")
                    tmp_atm_lines = natsorted(tmp_atm_lines)

                    for atmcnt, sorted_line in enumerate(tmp_atm_lines):
                        # atmcnt starts with 0
                        lmpdat_stuff, csitnam, cres = sorted_line
                        # translate original atom-ids to new internal ids
                        atm_id_old_new[int(lmpdat_stuff[0])] = atmcnt
                        #self.atm_idx_id[atmcnt] = int(lmpdat_stuff[0])

                        # check if instance of Atom with id atmcnt already exists
                        # i.e. a data file have had already been loaded
                        try:
                            cur_atm = self.atoms[atmcnt]

                            # overwrite data
                            if overwrite_data is True:
                                #cur_atm.atm_id  = int(lmpdat_stuff[0])
                                cur_atm.atm_id  = atmcnt
                                cur_atm.grp_id  = int(lmpdat_stuff[1])
                                cur_atm.atm_key = atm_tp_old_new[int(lmpdat_stuff[2])]
                                cur_atm.chge    = float(lmpdat_stuff[3])

                                # parse cgcmm stuff if available
                                if csitnam is not None:
                                    cur_atm.sitnam = csitnam

                                if cres is not None:
                                    cur_atm.res = cres

                            else:  # complement data

                                if not hasattr(self.atoms[atmcnt], "atm_id"):
                                    #cur_atm.atm_id  = int(lmpdat_stuff[0])
                                    cur_atm.atm_id  = atmcnt

                                if not hasattr(self.atoms[atmcnt], "grp_id"):
                                    cur_atm.grp_id  = int(lmpdat_stuff[1])

                                if not hasattr(self.atoms[atmcnt], "atm_key"):
                                    cur_atm.atm_key = atm_tp_old_new[int(lmpdat_stuff[2])]

                                if not hasattr(self.atoms[atmcnt], "chge"):
                                    cur_atm.chge    = float(lmpdat_stuff[3])

                                # parse cgcmm stuff if available
                                if not hasattr(self.atoms[atmcnt], "sitnam") and csitnam is not None:
                                    cur_atm.sitnam = csitnam

                                if not hasattr(self.atoms[atmcnt], "cres") and cres is not None:
                                    cur_atm.res = cres

                        # new atom must be created
                        except IndexError:
                            #atm_id=int(lmpdat_stuff[0])
                            cur_atm = mds.Atom(atm_id=atmcnt,
                                               grp_id=int(lmpdat_stuff[1]),
                                               atm_key=atm_tp_old_new[int(lmpdat_stuff[2])],
                                               chge=float(lmpdat_stuff[3])
                                               )

                            # # parse cgcmm stuff if available
                            if csitnam is not None:
                                cur_atm.sitnam = csitnam

                            if cres is not None:
                                cur_atm.res = cres

                            # append new atom if none was present before
                            self.atoms.append(cur_atm)

                        # parse coordinates
                        ccoords = np.array([float(i) for i in lmpdat_stuff[4:7]])
                        # append coordinates to temporary frame
                        tmp_ts.append(ccoords)

                    # append coordinates from data to (given) timesteps
                    self.ts_coords.append(np.array(tmp_ts))

                # /// bonds entry ///
                elif line.startswith("Bonds"):
                    lmpdat_in.next()  # skip empty line

                    # parse bonds
                    for bndcnt in xrange(total_bnds):
                        line = lmpdat_in.next()
                        line, comment = self._split_line(line)
                        line = line.split()
                        #line = [int(i) for i in line]

                        cur_bnd = mds.Bond()
                        #cur_bnd.bnd_id = int(line[0])
                        cur_bnd.bnd_id  = bndcnt
                        cur_bnd.bnd_key = bnd_tp_old_new[int(line[1])]
                        # translate original atom-ids
                        cur_bnd.atm_id1 = atm_id_old_new[int(line[2])]
                        cur_bnd.atm_id2 = atm_id_old_new[int(line[3])]

                        if comment is not None:
                            comment = comment.split()

                            # try reading the bond order if first item after
                            # the comment is a number
                            try:
                                cur_bnd.bnd_order = float(comment[0])

                                # read atom types the bond is between
                                if len(comment) > 2:
                                    cur_bnd.sitnam_2 = comment[2]
                                if len(comment) > 1:
                                    cur_bnd.sitnam_1 = comment[1]
                            except ValueError:
                                pass

                        self.bonds.append(cur_bnd)

                # /// angles entry ///
                elif "Angles" in line:
                    lmpdat_in.next()  # skip empty line

                    # parse angles
                    for angcnt in xrange(total_angs):
                        line = lmpdat_in.next()
                        line, comment = self._split_line(line)
                        line = line.split()
                        #line = [int(i) for i in line]

                        cur_ang = mds.Angle()
                        #cur_ang.ang_id = int(line[0])
                        cur_ang.ang_id  = angcnt
                        cur_ang.ang_key = ang_tp_old_new[int(line[1])]
                        cur_ang.atm_id1 = atm_id_old_new[int(line[2])]
                        cur_ang.atm_id2 = atm_id_old_new[int(line[3])]
                        cur_ang.atm_id3 = atm_id_old_new[int(line[4])]

                        self.angles.append(cur_ang)

                # /// dihedrals entry ///
                elif "Dihedrals" in line:
                    lmpdat_in.next()  # skip empty line

                    # parse dihedrals
                    for dihcnt in xrange(total_dihs):
                        line = lmpdat_in.next()
                        line, comment = self._split_line(line)
                        line = line.split()
                        #line = [int(i) for i in line]

                        cur_dih = mds.Dihedral()
                        #cur_dih.dih_id = int(line[0])
                        cur_dih.dih_id  = dihcnt
                        cur_dih.dih_key = dih_tp_old_new[int(line[1])]
                        cur_dih.atm_id1 = atm_id_old_new[int(line[2])]
                        cur_dih.atm_id2 = atm_id_old_new[int(line[3])]
                        cur_dih.atm_id3 = atm_id_old_new[int(line[4])]
                        cur_dih.atm_id4 = atm_id_old_new[int(line[5])]

                        self.dihedrals.append(cur_dih)

                # /// impropers entry ///
                elif "Impropers" in line:
                    lmpdat_in.next()  # skip empty line

                    for impcnt in xrange(total_imps):
                        line = lmpdat_in.next()
                        line, comment = self._split_line(line)
                        line = line.split()
                        #line = [int(i) for i in line]

                        cur_imp = mds.Improper()
                        #cur_imp.imp_id = int(line[0])
                        cur_imp.imp_id  = impcnt
                        cur_imp.imp_key = imp_tp_old_new[int(line[1])]
                        cur_imp.atm_id1 = atm_id_old_new[int(line[2])]
                        cur_imp.atm_id2 = atm_id_old_new[int(line[3])]
                        cur_imp.atm_id3 = atm_id_old_new[int(line[4])]
                        cur_imp.atm_id4 = atm_id_old_new[int(line[5])]

                        self.impropers.append(cur_imp)

                elif "Velocities" in line:
                    pass  # wip
                elif "# Forces" in line:
                    # since lammps also provides us with forces information,
                    # we make it possible to read those from the data file
                    # THIS IS NOT PART OF THE OFFICIAL LAMMPS DATA STRUCTURE!
                    lmpdat_in.next()  # skip empty line
                    tmp_forces = []

                    for atmcnt, sorted_line in enumerate(tmp_atm_lines):
                        line = lmpdat_in.next()
                        # parse coordinates
                        cforces = np.array([float(i) for i in line.split()[2:4]])
                        # append coordinates to temporary frame
                        tmp_forces.append(cforces)

                    # append coordinates from data to (given) timesteps
                    self.ts_forces.append(tmp_forces)
                else:
                    pass

        # append data-box to timestep-box (every timestep needs a box!)
        self.ts_boxes.append(lmpdat_box)

        # assign all atoms to their corresponding molecules internally
        self.fetch_molecules_by_bonds()

        # assign atoms to groups according to their molecule membership
        self.mols_to_grps()

        # only mix, if ii-pairs given (ij-pairs may be defined differently)
        if pair_ii is True:
            self.mix_pair_types(mode="ii", mix_style="arithmetic", debug=False)

        # more info
        if debug is True:
            print("***Info: Wiped all indices! Internal indices now start with 0.")

        # check charge of the system
        total_charge = sum([float(i.chge) for i in self.atoms])

        if debug is True:
            print("***Info: Total charge of the system is {}".format(total_charge))

    def _parse_cgcmm(self, cur_line):
        """
        Helper function for e.g. read_lmpdat to parse additional cgcmm-info
        from current line.
        """
        cgcmm_info  = None
        lmpdat_info = None
        sitnam      = None
        residue     = None

        # check if a comment is in line -> cgcmm-info (should be) in line
        if "#" in cur_line:
            line = cur_line.split("#")  # split by '#'
            lmpdat_info = line[0].split()
            cgcmm_info  = line[1].split()

            if len(cgcmm_info) > 0:
                sitnam  = cgcmm_info[0]

            if len(cgcmm_info) >= 2:
                residue = cgcmm_info[1]

        else:
            lmpdat_info = cur_line.split()

        return (lmpdat_info, sitnam, residue)

    def _split_line(self, curline):
        """
        Split line into comment and data.
        """
        data, comment = None, None

        if '#' in curline:
            line = curline.split("#")
            data = line[0]
            comment = line[1]
        else:
            data = curline

        return (data, comment)

    def write_lmpdat(self, lmpdat, frame_id=None, title=False, cgcmm=False):
        """
        Write new lmpdat.
        Sources:    http://lammps.sandia.gov/doc/2001/data_format.html
                    http://lammps.sandia.gov/doc/improper_style.html
        Further reading on string formatting:   https://pyformat.info/

        title:      str; title line of data file
        """
        # convert all internal indices back to original
        #for ccontainer in ("atms", "bnds", "angs", "dihs", "imps"):
        #    self._atm_bnd_ang_dih_imp_id_idx(entry=ccontainer, vice_versa=True)

        if title is False:
            title = "TITLE-LINE"

        # check box-type and convert if necessary
        if frame_id is None:
            frame_id = 0

        # convert (the copy) of the current box-type to a fractional box-type
        try:
            cbox = copy.copy(self.ts_boxes[frame_id])

            if cbox.boxtype == "cartesian":
                cbox.box_cart2lmp()
            elif cbox.boxtype == "lattice":
                cbox.box_lat2lmp()
            else:
                pass  # already of lammps' box-type
        except IndexError:
            cbox = None
            print("***Warning: No box specified - a simulation box must be specified for this file to work!")

        total_atms     = len(self.atoms)
        total_bnds     = len(self.bonds)
        total_angs     = len(self.angles)
        total_dihs     = len(self.dihedrals)
        total_imps     = len(self.impropers)
        total_atmtypes = len(self.atm_types)
        total_bndtypes = len(self.bnd_types)
        total_angtypes = len(self.ang_types)
        total_dihtypes = len(self.dih_types)
        total_imptypes = len(self.imp_types)
        longest_atm_id = len(str(len(self.atoms)))

        with open(lmpdat, "w") as lmpdat_out:

            # use cgcmm-style if wanted
            if cgcmm:
                lmpdat_out.write("{}; CGCMM style\n".format(title))
            else:
                lmpdat_out.write("{}\n".format(title))

            lmpdat_out.write("\n")
            # /// general summary ///
            lmpdat_out.write("{:>8d} atoms\n".format(total_atms))
            lmpdat_out.write("{:>8d} bonds\n".format(total_bnds))
            lmpdat_out.write("{:>8d} angles\n".format(total_angs))
            lmpdat_out.write("{:>8d} dihedrals\n".format(total_dihs))
            lmpdat_out.write("{:>8d} impropers\n".format(total_imps))
            lmpdat_out.write("\n")

            # /// force field summary ///
            lmpdat_out.write("{:>8d} atom types\n".format(total_atmtypes))
            lmpdat_out.write("{:>8d} bond types\n".format(total_bndtypes))
            lmpdat_out.write("{:>8d} angle types\n".format(total_angtypes))
            lmpdat_out.write("{:>8d} dihedral types\n".format(total_dihtypes))
            lmpdat_out.write("{:>8d} improper types\n".format(total_imptypes))
            lmpdat_out.write("\n")

            # /// box summary ///
            # write only if attributes were given beforehand

            # skip if no box was defined in data file
            if cbox is not None:
                if cbox.lmp_xlo is not None and cbox.lmp_xhi is not None:
                    lmpdat_out.write("{:> 12.6f} {:> 12.6f}  xlo xhi\n".format(
                        cbox.lmp_xlo, cbox.lmp_xhi)
                    )
                if cbox.lmp_ylo is not None and cbox.lmp_yhi is not None:
                    lmpdat_out.write("{:> 12.6f} {:> 12.6f}  ylo yhi\n".format(
                        cbox.lmp_ylo, cbox.lmp_yhi)
                    )
                if cbox.lmp_zlo is not None and cbox.lmp_zhi is not None:
                    lmpdat_out.write("{:> 12.6f} {:> 12.6f}  zlo zhi\n".format(
                        cbox.lmp_zlo, cbox.lmp_zhi)
                    )
                if cbox.lmp_xy is not None and cbox.lmp_xz is not None and cbox.lmp_yz is not None:
                    lmpdat_out.write("{:> 12.6f} {:> 12.6f} {:> 12.6f} xy xz yz\n".format(
                        cbox.lmp_xy, cbox.lmp_xz, cbox.lmp_yz)
                    )
                lmpdat_out.write("\n")

            # /// write atom types (masses) ///
            if self.atm_types:
                lmpdat_out.write("Masses\n")
                lmpdat_out.write("\n")

                for iatyp in sorted(self.atm_types):
                    lmpdat_out.write("{:>8d} {:>12.4f} ".format(
                        iatyp, self.atm_types[iatyp].weigh)
                    )

                    # write atom-name if existent
                    if cgcmm:
                        try:
                            lmpdat_out.write("# {}".format(self.atm_types[iatyp].sitnam))
                        except AttributeError:
                            pass

                    lmpdat_out.write("\n")
                lmpdat_out.write("\n")

            # /// write bond types(coeffs) ///
            if self.bnd_types:
                lmpdat_out.write("Bond Coeffs\n")
                lmpdat_out.write("\n")

                for ibtyp in self.bnd_types:
                    lmpdat_out.write("{:>8d} {:>12.6f} {:>12.6f}".format(
                        ibtyp,
                        self.bnd_types[ibtyp].prm1,
                        self.bnd_types[ibtyp].prm2)
                    )

                    # write further parameters if existent
                    try:
                        if self.bnd_types[ibtyp].prm3:
                            lmpdat_out.write("{:>12.6f}".format(
                                self.bnd_types[ibtyp].prm3)
                            )
                        if self.bnd_types[ibtyp].prm4:
                            lmpdat_out.write("{:>12.6f}".format(
                                self.bnd_types[ibtyp].prm4)
                            )
                    except AttributeError:
                        pass

                    if self.bnd_types[ibtyp].comment is not None:
                        lmpdat_out.write("  #{}".format(self.bnd_types[ibtyp].comment))

                    lmpdat_out.write("\n")
                lmpdat_out.write("\n")

            # /// write angle types(coeffs) ///
            if self.ang_types:
                lmpdat_out.write("Angle Coeffs\n")
                lmpdat_out.write("\n")

                for iangtyp in self.ang_types:
                    lmpdat_out.write("{:>8d} {:>12.6f} {:>12.6f}".format(
                        iangtyp, self.ang_types[iangtyp].prm1, self.ang_types[iangtyp].prm2)
                    )

                    # write further parameters if existent
                    try:
                        if self.ang_types[iangtyp].prm3:
                            lmpdat_out.write("{:>12.6f}".format(
                                self.ang_types[iangtyp].prm3)
                            )
                        if self.ang_types[iangtyp].prm4:
                            lmpdat_out.write("{:>12.6f}".format(
                                self.ang_types[iangtyp].prm4)
                            )
                    except AttributeError:
                        pass

                    if self.ang_types[iangtyp].comment is not None:
                        lmpdat_out.write("  #{}".format(self.ang_types[iangtyp].comment))

                    lmpdat_out.write("\n")
                lmpdat_out.write("\n")

            # /// dihedral types(coeffs) ///
            if self.dih_types:
                lmpdat_out.write("Dihedral Coeffs\n")
                lmpdat_out.write("\n")

                for idtyp in self.dih_types:
                    try:
                        if not self.dih_types[idtyp].weigh_factor:
                            self.dih_types[idtyp].weigh_factor = 0
                    except AttributeError:
                        self.dih_types[idtyp].weigh_factor = 0

                    lmpdat_out.write("{:>8d} {:>12.6f} {:>5d} {:>6d} {:>15.6f}".format(
                        idtyp,
                        self.dih_types[idtyp].prm_k,
                        int(self.dih_types[idtyp].prm_n),
                        int(self.dih_types[idtyp].prm_d),
                        self.dih_types[idtyp].weigh_factor)
                    )

                    if self.dih_types[idtyp].comment is not None:
                        lmpdat_out.write("  #{}".format(self.dih_types[idtyp].comment))

                    lmpdat_out.write("\n")
                lmpdat_out.write("\n")

            # /// improper types(coeffs) ///
            if self.imp_types:
                lmpdat_out.write("Improper Coeffs\n")
                lmpdat_out.write("\n")

                for iityp in self.imp_types:
                    lmpdat_out.write("{:>8d} {:>12.6f} {:>5d} {:>6d}".format(
                        iityp,
                        self.imp_types[iityp].prm_k,
                        int(self.imp_types[iityp].prm_d),
                        int(self.imp_types[iityp].prm_n))
                    )

                    if self.imp_types[iityp].comment is not None:
                        lmpdat_out.write("  #{}".format(self.imp_types[iityp].comment))

                    lmpdat_out.write("\n")
                lmpdat_out.write("\n")

            # /// pair types (coeffs) (requires mixed vdw)///
            if hasattr(self, "pair_types") is True:
                if self.pair_types != []:
                    # crawl through pair_types if there are any ij-pairs
                    pair_ij = False
                    for cpair in self.pair_types:
                        if cpair.pairs == "ij":
                            pair_ij = True
                            break

                    if pair_ij is True:
                        lmpdat_out.write("PairIJ Coeffs\n")
                    else:
                        lmpdat_out.write("Pair Coeffs\n")

                    lmpdat_out.write("\n")

                    for prtyp in self.pair_types:

                        if prtyp.pairs == "ii":
                            lmpdat_out.write("{:>8d} {:>12.6f} {:>12.6f}".format(
                                prtyp.atm_key_i,
                                prtyp.epsilon_ij,
                                prtyp.sigma_ij))
                        elif prtyp.pairs == "ij":
                            lmpdat_out.write("{:>8d} {:>8d} {:>12.6f} {:>12.6f}".format(
                                prtyp.atm_key_i,
                                prtyp.atm_key_j,
                                prtyp.epsilon_ij,
                                prtyp.sigma_ij))
                        else:
                            raise Warning("No type for current pair given!")

                        # add some comment stuff
                        try:
                            lmpdat_out.write(" # {}".format(prtyp.lr_key))
                        except AttributeError:
                            pass

                        lmpdat_out.write("\n")

                    lmpdat_out.write("\n")

            # /// atoms entry ///
            if self.atoms:
                lmpdat_out.write("Atoms\n")
                lmpdat_out.write("\n")
                try:
                    longest_grp_id = len(str(self.atoms[-1].grp_id))
                except AttributeError:
                    longest_grp_id = 2

                longest_atm_key = len(str(len(self.atm_types)))

                for cidx, catm in enumerate(self.atoms):
                    # get atom-index for current atom
                    #cidx = self.atm_id_idx[catm.atm_id]
                    #cidx = catm.atm_id
                    # atom-id
                    if not hasattr(catm, "atm_id"):
                        catm.atm_id = cidx

                    if not hasattr(catm, "grp_id"):
                        catm.grp_id = 1

                    if not hasattr(catm, "atm_key"):
                        catm.atm_key = 1

                    if not hasattr(catm, "chge"):
                        catm.chge = 0.0

                    lmpdat_out.write("{0:<8d} {1:<{width_2}d}      {2:<{width_3}d} {3: >10.6f} {c[0]: >16.6f} {c[1]: >12.6f} {c[2]: >12.6f}".format(
                        catm.atm_id,
                        catm.grp_id,
                        catm.atm_key,
                        catm.chge,
                        width_2=longest_grp_id,
                        width_3=longest_atm_key,
                        c=self.ts_coords[frame_id][cidx])
                    )

                    # write cgcmm info (if given)
                    if cgcmm:
                        lmpdat_out.write(" #")
                        try:
                            lmpdat_out.write(" {:<s}".format(catm.sitnam))
                            lmpdat_out.write(" {}".format(catm.res))
                        except AttributeError:
                            pass

                    lmpdat_out.write("\n")
                lmpdat_out.write("\n")

            # /// bonds entry ///
            if self.bonds:
                lmpdat_out.write("Bonds\n")
                lmpdat_out.write("\n")
                longest_bnd_key = len(str(len(self.bnd_types)))

                for cbnd in self.bonds:
                    lmpdat_out.write("{0:<8d} {1:>{width_2}d}      {2:>{width_3}d} {3:>{width_3}d}".format(
                        cbnd.bnd_id, cbnd.bnd_key, cbnd.atm_id1, cbnd.atm_id2,
                        width_2=longest_bnd_key,
                        width_3=longest_atm_id)
                    )

                    # write bond order as well if given for current bond
                    if hasattr(cbnd, "bnd_order"):
                        lmpdat_out.write(" # {}".format(cbnd.bnd_order))


                    lmpdat_out.write("\n")
                lmpdat_out.write("\n")

            # /// angles entry ///
            if self.angles:
                lmpdat_out.write("Angles\n")
                lmpdat_out.write("\n")
                longest_ang_key = len(str(len(self.ang_types)))

                for cang in self.angles:
                    lmpdat_out.write("{0:<8d} {1:>{width_2}d}      {2:>{width_3}d} {3:>{width_3}d} {4:>{width_3}d}".format(
                        cang.ang_id, cang.ang_key,
                        cang.atm_id1, cang.atm_id2, cang.atm_id3,
                        width_2=longest_ang_key,
                        width_3=longest_atm_id)
                    )

                    lmpdat_out.write("\n")
                lmpdat_out.write("\n")

            # /// dihedrals entry ///
            if self.dihedrals:
                lmpdat_out.write("Dihedrals\n")
                lmpdat_out.write("\n")
                longest_dih_key = len(str(len(self.dih_types)))

                for cdih in self.dihedrals:
                    lmpdat_out.write("{0:<8d} {1:>{width_2}d}      {2:>{width_3}d} {3:>{width_3}d} {4:>{width_3}d} {5:>{width_3}d}".format(
                        cdih.dih_id, cdih.dih_key,
                        cdih.atm_id1, cdih.atm_id2, cdih.atm_id3, cdih.atm_id4,
                        width_2=longest_dih_key,
                        width_3=longest_atm_id)
                    )
                    lmpdat_out.write("\n")
                lmpdat_out.write("\n")

            if self.impropers:
                lmpdat_out.write("Impropers\n")
                lmpdat_out.write("\n")
                longest_imp_key = len(str(len(self.imp_types)))

                for cimp in self.impropers:
                    lmpdat_out.write("{0:<8d} {1:>{width_2}d}      {2:>{width_3}d}  {3:>{width_3}d} {4:>{width_3}d} {5:>{width_3}d}".format(
                        cimp.imp_id,
                        cimp.imp_key,
                        cimp.atm_id1, cimp.atm_id2, cimp.atm_id3, cimp.atm_id4,
                        width_2=longest_imp_key,
                        width_3=longest_atm_id)
                    )
                    lmpdat_out.write("\n")
                lmpdat_out.write("\n")

            if self.ts_forces:
                lmpdat_out.write("# Forces\n\n")

                for cidx, catm in enumerate(self.atoms):
                    lmpdat_out.write("# {0:<8d} {c[0]: >16.6f} {c[1]: >12.6f} {c[2]: >12.6f}".format(
                        catm.atm_id,
                        c=self.ts_forces[frame_id][cidx])
                    )
                    lmpdat_out.write("\n")
                lmpdat_out.write("\n")

    def import_dcd(self, dcd):
        """
        Read and write the "DCD" binary trajectory file format used by LAMMPS.
        Other formats are not considered (yet).
        Sources:    http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
                    https://docs.python.org/2/library/struct.html#functions-and-exceptions
                    http://www.devdungeon.com/content/working-binary-data-python
                    http://stackoverflow.com/questions/1035340/reading-binary-file-in-python-and-looping-over-each-byte/1035360#1035360
                    http://stackoverflow.com/questions/38297929/why-does-unpacking-a-struct-result-in-a-tuple
                    http://prody.csb.pitt.edu/_modules/prody/trajectory/dcdfile.html#codemodal
        Open DCD and read the remarks. This function must be called before
        anything else that has to do anything with file reading.
        """
        self._dcdfile = open(dcd, 'rb')
        self._read_header()
        self._read_title()
        # get position in file (bytes?) after header and title
        self._pos_1 = self._dcdfile.tell()

    def jump_to_first_frame(self):
        """
        Rewind to the first frame
        """
        self._dcdfile.seek(self._pos_1)

    def _read_header(self, debug=False):
        """
        Read header block.
        """
        self.extra_blck = None  # only if unit cell
        self.has_4dims  = None  # purpose unknown
        self.is_charmm  = False

        # header block
        hdr_blck = agldh.read_record(self._dcdfile)
        hdr = struct.unpack('4c9if10i', hdr_blck)
        self.nframes = hdr[4]  # total number of frames
        self.sframe  = hdr[5]  # number of start frame
        self.step    = hdr[6]  # number of frames between each frame
        self.lframe  = hdr[7]  # number of last frame

        # check charmm-formatting
        if hdr[23] != 0:
            self.is_charmm = True
            self.extra_blck = hdr[14]
            # some dcd files have 4 dimensions?
            self.has_4dims = hdr[15]

            if debug is True:
                print("***Reading Charmm formatted DCD with 4 dimensions!***")

    def _read_title(self, debug=False):
        """
        Read title blocks (only possible if header block was read before!)
        """
        # 2 title lines
        title_1 = agldh.read_record(self._dcdfile)  # 1st title block
        title_2 = agldh.read_record(self._dcdfile)  # 2nd title block
        self.natoms, = struct.unpack("i", title_2)

        if debug is True:
            print("   Remark 1: {}\n   Remark 2: Number of Atoms: {}".format(
                  title_1, self.natoms))

    def _read_frame(self):
        """
        Read one frame.
        Layout of unitcell is [A, alpha, B, beta, gamma, C] (Historical reasons)
        """
        # read cell information
        if self.extra_blck:
            cur_blck = agldh.read_record(self._dcdfile)
            cur_cell = struct.unpack("6d", cur_blck)

        # read coordinate-sets
        x_coordset = agldh.read_record(self._dcdfile)
        y_coordset = agldh.read_record(self._dcdfile)
        z_coordset = agldh.read_record(self._dcdfile)
        x = np.fromstring(x_coordset, dtype=np.dtype('f'), count=self.natoms)
        y = np.fromstring(y_coordset, dtype=np.dtype('f'), count=self.natoms)
        z = np.fromstring(z_coordset, dtype=np.dtype('f'), count=self.natoms)

        # 4th dimension given? (has also to be read)
        if self.has_4dims:
            agldh.read_record(self._dcdfile)
            #dims_4_blck = agldh.read_record(self._dcdfile)
            #dunno = struct.unpack(str(self.natoms)+"i", dims_4_blck)

        return(x, y, z, cur_cell)

    def _skip_frame(self):
        if self.extra_blck:
            agldh.read_record(self._dcdfile)

        agldh.read_record(self._dcdfile)
        agldh.read_record(self._dcdfile)
        agldh.read_record(self._dcdfile)

        if self.has_4dims:
            agldh.read_record(self._dcdfile)

    def read_frames(self, frame=None, to_frame=-1, frame_by="index", debug=False):
        """
        Sources:    https://github.com/MDAnalysis/mdanalysis/issues/187
        """
        if debug:
            print("***Verbose: MB already read: {:.2f} MiB".format(self._dcdfile.tell() / 1000000))

        M_PI_2 = np.pi / 2
        # convert input to corresponding indices
        frm, to_frm = agldh.reshape_arguments(self.sframe, self.nframes,
                                              self.step, frame, to_frame,
                                              frame_by)

        # preallocate memory for arrays to come
        if to_frame is None:
            num_frames = abs(frm - to_frm) - 1
        else:
            num_frames = abs(frm - to_frm)  # one frame is always read

        if debug is True:
            print("***Info: Reading: Frame (start): {}, ToFrame (excluded): {}, NumFrames: {}".format(frm, to_frm, num_frames))

        x_set, y_set, z_set = agldh.deploy_array(num_frames, self.natoms)
        # fill arrays with coordinates
        ptr = 0  # pointer to place data in right position of array

        for frame_num in xrange(self.nframes):

            if frm <= frame_num < to_frm:

                # coords
                x, y, z, cur_box = self._read_frame()
                x_set[ptr] = x
                y_set[ptr] = y
                z_set[ptr] = z

                # create box and append to other boxes of trajectory
                #TODO: Check if angles are right this way with triclinic cell
                alpha = np.radians(90.0 - np.arcsin(cur_box[4])*90.0/M_PI_2)  # cosAB
                beta  = np.radians(90.0 - np.arcsin(cur_box[3])*90.0/M_PI_2)  # cosAC
                gamma = np.radians(90.0 - np.arcsin(cur_box[1])*90.0/M_PI_2)  # cosBC
                cur_box = mdb.Box(ltc_alpha=alpha,
                                  ltc_beta=beta,
                                  ltc_gamma=gamma,
                                  ltc_a=cur_box[0],
                                  ltc_b=cur_box[2],
                                  ltc_c=cur_box[5],
                                  boxtype="lattice")

                # convert lattice-box-type to lammps-box-type
                cur_box.box_lat2lmp()
                self.ts_boxes.append(cur_box)

                # get step numbers per frame so we can later access them if wanted
                #self.ts_steps.append(frame_num*self.step+self.sframe)
                ptr += 1
            elif frame_num > to_frm:
                break
            else:
                self._skip_frame()  # keep reading until first is reached

        # concatenate all arrays to (1,3)-arrays
        coordinates = np.stack((x_set, y_set, z_set), axis=-1)

        # append coordinates to universe ts-coordinates
        for i in coordinates:
            self.ts_coords.append(i)

    def append_dcds(self, *dcd_files, **read_frame_args):
        """
        Read and append another DCD to an existing one/s.
        Only dcd may appended, if the number of atoms is the same!

        dcds:   str; dcd file name(s) to append
        read_frame_args: same as for method read_frames
        """
        cur_num_atoms = self.natoms  # same atoms as in data

        for cur_dcd in dcd_files:
            c_dcd = LmpStuff()
            c_dcd.import_dcd(cur_dcd)

            # check if number of atoms is the same as in data/previous dcd file
            if cur_num_atoms != c_dcd.natoms:
                raise RuntimeError("Different number of atoms in DCD-files!")

            cur_num_atoms = c_dcd.natoms
            # use same kwargs as for "append_dcds"
            c_dcd.read_frames(**read_frame_args)
            do_not_append = []

            # find duplicate step-entries, save indices
            for k, i in enumerate(c_dcd.ts_coords):
                for j in self.ts_coords:
                    if i == j:
                        do_not_append.append(k)

            # append ts-steps, ts-boxes and ts-coordinates to universe
            for iidx, istp in enumerate(c_dcd.ts_coords):
                if iidx not in do_not_append:
                    self.ts_coords.append(istp)
                    self.ts_boxes.append(c_dcd.ts_boxes[iidx])
                    self.ts_coords.append(c_dcd.ts_coords[iidx])

            c_dcd.close_dcd()

    def close_dcd(self, debug=False):
        """
        Close dcd-file if still open.
        """
        if not self._dcdfile.closed:

            if debug is True:
                print("***Info: Closing file: {}.".format(self._dcdfile))

            self._dcdfile.close()

    def write_dcd(self, *frames):
        """
        Write frames to DCD-file. TBD.
        """
        pass


################################################################################
# Shortcut functions for common procedures
################################################################################

def read_lmpdat(lmpdat, dcd=None, frame_idx_start=-1, frame_idx_stop=-1):
    """
    Read a lammps data file and optionally a dcd file on top.

    This function simplifies the process of reading a lammps data file and loading
    a dcd file on top of it for further coordinates. It returns a Universe object
    which has many methods to manipulate it with.

    Parameters
    ----------
    lmpdat : str
        Name of the lammps data file

    dcd : str (optional)
        Name of the dcd file with the same amounts of atoms as the lammps data
        file

    frame_idx : int (default: -1)
        Index of the frame to use from the dcd file.

    Returns
    -------
    md_sys : LmpStuff object
        An object of LmpStuff which can be further processed.

    """
    # read output from quenching
    md_sys = LmpStuff()
    md_sys.read_lmpdat(lmpdat)

    if dcd is not None:
        md_sys.import_dcd(dcd)
        # since we are only interested in one frame, delete all others
        md_sys.ts_coords = []
        md_sys.ts_boxes = []

        # enable reading the last frame with negative indexing
        if frame_idx_start == frame_idx_stop:
            if frame_idx_stop < 0:
                md_sys.read_frames(frame=frame_idx_stop - 1, to_frame=frame_idx_stop)
            else:
                md_sys.read_frames(frame=frame_idx_stop, to_frame=frame_idx_stop + 1)
        else:
            if frame_idx_stop == -1:
                md_sys.read_frames(frame=frame_idx_start, to_frame=frame_idx_stop)
            else:
                md_sys.read_frames(frame=frame_idx_start, to_frame=frame_idx_stop + 1)

        md_sys.close_dcd()

    return md_sys


def write_lmpdat(lmpdat_out, lmpdat_a, lmpdat_b=None, dcd_a=None, dcd_b=None, frame_idx_a=-1, frame_idx_b=-1, pair_coeffs=None):
    """
    Write a lammps data file.

    Write a lammps data file by reading the file and merging it with a second
    lammps data file (optional). DCD files may be loaded as well for each system
    on top. When merging, simulation boxes from system a will be overwritten by
    the simulation boxes from system b.

    Parameters
    ----------
    lmpdat_out : str
        name of the lammps output file

    lmpdat_a : str
        lammps data file a

    lmpdat_b : str (optional)
        lammps data file b

    dcd_a : str
        dcd file which has the same amount of coordinates as lammps data a has atoms


    dcd_b : str (optional)
        dcd file which has the same amount of coordinates as lammps data b has atoms

    frame_idx_a : int (default: -1)
        index of the frame to use from dcd a (default: last frame)

    frame_idx_b : int (default : -1)
        index of the frame to use from dcd b (default: last frame)

    pair_coeffs : str (optional)
        mixing type for pair coefficients (default: None): 'ii' or 'jj'

    """
    sys_a = read_lmpdat(lmpdat_a, dcd_a, frame_idx_a, frame_idx_a)

    if lmpdat_b is not None:
        sys_b = read_lmpdat(lmpdat_b, dcd_b, frame_idx_b, frame_idx_b)
        # since we only read one frame, only this frame combined will be written
        sys_ab = mdu.merge_systems([sys_a, sys_b], pair_coeffs=pair_coeffs)
    else:
        sys_ab = sys_a

    sys_ab.change_indices(incr=1, mode="increase")
    sys_ab.write_lmpdat(lmpdat_out, cgcmm=True)


def cut_box(lmpdat_out, lmpdat, box, dcd=None, frame_idx=-1):
    """
    Cut a smaller solvent box from a bigger one that will be used during the simulation.

    Given a (in the best case) larger solvent box, a smaller one will be cut.
    Having only as many solvent molecules as absolutely necessary reduces the
    calculation time. This is possible since the potential energy for a group
    of atoms may be calculated with lammps. The center of the box will be at
    the origin.
    Only reads one frame.

    Parameters
    ----------
    lmpdat_out : str
        name of new lammps data file with cut coordinates

    lmpdat : str
        lammps data file of the system to cut from

    box : Box
        box to cut out from lmpdat

    dcd : str (optional)
        Name of the dcd file with the same amounts of atoms as the lammps data
        file

    frame_idx : int (default: -1)
        Index of the frame to use from the dcd file.

    """
    md_sys = read_lmpdat(lmpdat, dcd, frame_idx_start=frame_idx - 1, frame_idx_stop=frame_idx)

    # generate planes that describe the box
    cart_box = copy.deepcopy(box)
    cart_box.box_lmp2cart()

    # define planes which form the box to cut (plane vectors always need the same
    # origin or they will be shifted)
    plane_ab = agv.get_plane(cart_box.crt_a, cart_box.crt_b)
    plane_ca = agv.get_plane(cart_box.crt_c, cart_box.crt_a)
    plane_bc = agv.get_plane(cart_box.crt_b, cart_box.crt_c)

    # opposite planes to ab, ca and bc
    plane_ab_c = [-1 * i for i in agv.get_plane(cart_box.crt_a, cart_box.crt_b, cart_box.crt_c)]
    plane_ca_b = [-1 * i for i in agv.get_plane(cart_box.crt_c, cart_box.crt_a, cart_box.crt_b)]
    plane_bc_a = [-1 * i for i in agv.get_plane(cart_box.crt_b, cart_box.crt_c, cart_box.crt_a)]

    # cut plane
    md_sys.cut_shape(-1, True, plane_ab, plane_ab_c, plane_ca, plane_ca_b, plane_bc, plane_bc_a)

    # buffer for box lengths to prevent clashes using pbc
    cart_box.box_cart2lat()
    cart_box.ltc_a += 4
    cart_box.ltc_b += 4
    cart_box.ltc_c += 4
    cart_box.box_lat2lmp()

    # def new box vectors by given box angles
    md_sys.ts_boxes[0] = cart_box

    if box.lmp_xy < 1e-10 and box.lmp_xz < 1e-10 and box.lmp_yz < 1e-10:
        md_sys.ts_boxes[0].lmp_xy = None
        md_sys.ts_boxes[0].lmp_xz = None
        md_sys.ts_boxes[0].lmp_yz = None

    md_sys.mols_to_grps()
    md_sys.change_indices(incr=1, mode="increase")
    md_sys.write_lmpdat(lmpdat_out, cgcmm=True)


#class LmpDihedral(object):
#    """Calculating the energy of the dihedral potential as in the lammps manual.
#    Source: dihedral_style charmm.
#    """
#    def __init__(self, k, n, d):
#        """
#        Set parameters.
#
#        K (energy)
#        n (multiplicity, >= 0)
#        d (degrees) = phi0 + 180
#
#        """
#        self.k = k
#        self.n = n
#        self.d = d
#        self.dihedral_energy = None
#
#    def calc(self, phi):
#        """
#        Calculate the value for phi given certain parameters.
#
#        phi : float
#            angle phi in degrees
#        """
#        self.dihedral_energy = self.k * (1 + np.cos(self.n * np.radians(phi) - np.radians(self.d)))
