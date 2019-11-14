
import math
import numpy as np
import md_stars as mds
import md_universe as mdu
import ag_prmtop_helper_functions as agphf

__version__ = "2017-08-17"


class AmberStuff(mdu.Universe):
    """
    Sources:    http://ambermd.org/formats.html
    """
    def __init__(self):
        """
        Get our containers ready.
        """
        mdu.Universe.__init__(self)

    def read_prmtop(self, prmtop):
        """
        Read the contents of the amber prmtop-file. CHARMM-Entries will not be
        read!

        Input:
            >   mode        str; complement|overwrite|append;
                            complement (missing) data, overwrite given data
                            or append to given data
        """
        # /// parse file
        with open(prmtop, "r") as prmtop_in:
            prmtop_version = prmtop_in.readline()

            # parse sections
            for line in prmtop_in:

                if line.startswith("%FLAG TITLE"):
                    # section contains the title of the topology file
                    next(prmtop_in)  # line with formatting info
                    prmtop_title = next(prmtop_in)
                elif line.startswith("%FLAG POINTERS"):
                    # section which contains the information about how many
                    # parameters are present in all of the sections

                    next(prmtop_in)  # line with formatting info

                    line = prmtop_in.next().split()
                    (natom, ntypes, nbonh, mbona, ntheth, mtheta, nphih, mphia,
                     nhparm, nparm) = [int(i) for i in line]

                    line = prmtop_in.next().split()
                    (nnb, nres, nbona, ntheta, nphia, numbnd, numang, nptra,
                     natyp, nphb) = [int(i) for i in line]

                    line = prmtop_in.next().split()
                    (ifpert, nbper, ngper, ndper, mbper, mgper, mdper, ifbox,
                     nmxrs, ifcap) = [int(i) for i in line]

                    line = prmtop_in.next().split()
                    numextra = line[0]

                    try:
                        ncopy = line[1]  # entry "copy" need not to be given
                    except(IndexError):
                        pass

                elif line.startswith("%FLAG ATOM_NAME"):
                    # SECTION: ATOM_NAME
                    # section which contains the atom name for every atom in
                    # the prmtop
                    section_atom_name = agphf.parse_section(prmtop_in, natom)

                elif line.startswith("%FLAG CHARGE"):
                    # SECTION: CHARGE
                    # section which contains the charge for every atom in the prmtop
                    # prmtop-charges must be divided by 18.2223 (reasons unknown)
                    section_charge = agphf.parse_section(prmtop_in, natom)

                elif line.startswith("%FLAG ATOMIC NUMBER"):
                    pass
                elif line.startswith("%FLAG MASS"):
                    # SECTION: MASS
                    # section which  contains the atomic mass of every atom
                    # in g/mol.
                    section_mass = agphf.parse_section(prmtop_in, natom)

                elif line.startswith("%FLAG ATOM_TYPE_INDEX"):
                    # SECTION ATOM TYPE INDEX
                    # section which contains the Lennard-Jones atom type index.
                    # The Lennard-Jones potential contains parameters for every
                    # pair of atoms in the system. All atoms with the same
                    # sigma and epsilon parameters are assigned to the same type (regardless
                    # of whether they have the same AMBER ATOM TYPE).
                    section_atom_type_index = agphf.parse_section(
                        prmtop_in, natom, itype="int")

                elif line.startswith("%FLAG NUMBER_EXCLUDED_ATOMS"):
                    # section which contains the number of atoms that need to be
                    # excluded from the non-bonded calculation loop for atom i
                    # because i is involved in a bond, angle, or torsion with those atoms
                    section_number_excluded_atoms = agphf.parse_section(
                        prmtop_in, natom, itype="int")

                elif line.startswith("%FLAG NONBONDED_PARM_INDEX"):
                    # section which contains the pointers for each pair of LJ
                    # atom types into the LENNARD JONES ACOEF and
                    # LENNARD JONES BCOEF arrays
                    section_nonbonded_parm_index = agphf.parse_section(
                        prmtop_in, ntypes*ntypes, itype="int")

                elif line.startswith("%FLAG RESIDUE_LABEL"):
                    # section which contains the residue name for every residue
                    # in the prmtop
                    section_residue_label = agphf.parse_section(
                        prmtop_in, nres)
                elif line.startswith("%FLAG RESIDUE_POINTER"):
                    section_residue_pointer = agphf.parse_section(
                        prmtop_in, nres, itype="int")
                elif line.startswith("%FLAG BOND_FORCE_CONSTANT"):
                    # section which lists all of the bond force constants
                    # k in kcal/(mol*Angstrom**2) for each unique bond type
                    section_bond_force_constant = agphf.parse_section(
                        prmtop_in, numbnd)
                elif line.startswith("%FLAG BOND_EQUIL_VALUE"):
                    # section which lists all of the bond equilibrium distances
                    # in "Angstrom" for each unique bond type
                    section_bond_equil_value = agphf.parse_section(
                        prmtop_in, numbnd)
                elif line.startswith("%FLAG ANGLE_FORCE_CONSTANT"):
                    # section which contains all of the angle equilibrium angles
                    # in "radians"
                    section_angle_force_constant = agphf.parse_section(
                        prmtop_in, numang)
                elif line.startswith("%FLAG ANGLE_EQUIL_VALUE"):
                    # section which contains all of the angle equilibrium angles
                    # in "radians"
                    section_angle_equil_value = agphf.parse_section(
                        prmtop_in, numang)
                elif line.startswith("%FLAG DIHEDRAL_FORCE_CONSTANT"):
                    # section which lists the torsion force constants in kcal/mol
                    # for each unique torsion type
                    section_dihedral_force_constant = agphf.parse_section(
                        prmtop_in, nptra)
                elif line.startswith("%FLAG DIHEDRAL_PERIODICITY"):
                    # section which lists the periodicity n for each unique
                    # torsion type; only int
                    section_dihedral_periodicity = agphf.parse_section(
                        prmtop_in, nptra)
                elif line.startswith("%FLAG DIHEDRAL_PHASE"):
                    # section which lists the phase shift for each unique
                    # torsion type in "radians"
                    section_dihedral_phase = agphf.parse_section(
                        prmtop_in, nptra)
                elif line.startswith("%FLAG SCEE_SCALE_FACTOR"):
                    # section which lists the factor by which 1-4 electrostatic
                    # interactions are divided (i.e., the two atoms on either
                    # end of a torsion)
                    section_scee_scale_factor = agphf.parse_section(
                        prmtop_in, nptra)
                elif line.startswith("%FLAG SCNB_SCALE_FACTOR"):
                    # section which lists the factor by which 1-4 van der Waals
                    # interactions are divided (i.e., the two atoms on either
                    # end of a torsion)
                    section_scnb_scale_factor = agphf.parse_section(
                        prmtop_in, nptra)
                elif line.startswith("%FLAG SOLTY"):
                    # section which is unused
                    section_solty = agphf.parse_section(prmtop_in, natyp)
                elif line.startswith("%FLAG LENNARD_JONES_ACOEF"):
                    # section contains the LJ A-coefficients for all pairs of
                    # distinct LJ types
                    section_lennard_jones_acoef = agphf.parse_section(
                        prmtop_in, ntypes*(ntypes+1)/2)
                elif line.startswith("%FLAG LENNARD_JONES_BCOEF"):
                    # section contains the LJ A-coefficients for all pairs of
                    # distinct LJ types
                    section_lennard_jones_bcoef = agphf.parse_section(
                        prmtop_in, ntypes*(ntypes+1)/2)
                elif line.startswith("%FLAG BONDS_INC_HYDROGEN"):
                    # section which contains a list of every bond in the system
                    # in which at least one atom is Hydrogen
                    section_bonds_inc_hydrogen = agphf.parse_section(
                        prmtop_in, 3*nbonh, chunksize=3, itype="int")
                elif line.startswith("%FLAG BONDS_WITHOUT_HYDROGEN"):
                    # section contains a list of every bond in the system in
                    # which neither atom is a Hydrogen
                    section_bonds_without_hydrogen = agphf.parse_section(
                        prmtop_in, 3*nbona, chunksize=3, itype="int")
                elif line.startswith("%FLAG ANGLES_INC_HYDROGEN"):
                    # section contains a list of every angle in the system in
                    # which at least one atom is Hydrogen
                    section_angles_inc_hydrogen = agphf.parse_section(
                        prmtop_in, 4*ntheth, chunksize=4, itype="int")
                elif line.startswith("%FLAG ANGLES_WITHOUT_HYDROGEN"):
                    # section which contains a list of every angle in the system
                    # in which no atom is Hydrogen
                    section_angles_without_hydrogen = agphf.parse_section(
                        prmtop_in, 4*ntheta, chunksize=4, itype="int")
                elif line.startswith("%FLAG DIHEDRALS_INC_HYDROGEN"):
                    # section contains a list of every torsion in the system in
                    # which at least one atom is Hydrogen
                    section_dihedrals_inc_hydrogen = agphf.parse_section(
                        prmtop_in, 5*nphih, chunksize=5, itype="int")
                elif line.startswith("%FLAG DIHEDRALS_WITHOUT_HYDROGEN"):
                    section_dihedrals_without_hydrogen = agphf.parse_section(
                        prmtop_in, 5*nphia, chunksize=5, itype="int")
                elif line.startswith("%FLAG AMBER_ATOM_TYPE"):
                    section_amber_atom_type = agphf.parse_section(
                        prmtop_in, natom)
                elif line.startswith("%FLAG TREE_CHAIN_CLASSIFICATION"):
                    agphf.entry_not_read("TREE_CHAIN_CLASSIFICATION")
                elif line.startswith("%FLAG JOIN_ARRAY"):
                    agphf.entry_not_read("JOIN_ARRAY")
                elif line.startswith("%FLAG IROTAT"):
                    agphf.entry_not_read("IROTAT")
                elif line.startswith("%FLAG RADIUS_SET"):
                    agphf.entry_not_read("RADIUS_SET")
                elif line.startswith("%FLAG RADII"):
                    section_radii = agphf.parse_section(
                        prmtop_in, natom)
                elif line.startswith("%FLAG SCREEN"):
                    agphf.entry_not_read("SCREEN")
                else:
                    pass

        # /// arrange data - force field section ///
        # /// atom-types
        #     (atoms with same epsilon and sigma (lj) have same type)

        # 1. get acoef/bcoef-ids for ii-interactions
        abcoef_ids = agphf.get_AB_ids(ntypes, section_nonbonded_parm_index)

        # 2. calculate sigma from coefs
        # dict to translate between internal and amber indices
        atm_key_old_new = {}

        for cidx, cid in enumerate(abcoef_ids):
            #cidx += 1
            cur_sig, cur_eps = agphf.sig_eps_from_AB(section_lennard_jones_acoef[cid],
                                                     section_lennard_jones_bcoef[cid])

            # skip duplicates
            try:
                if self.atm_types[cidx]:
                    print("***Prmtop-Info: Skipping atom type {}".format(cid))
            except KeyError:
                self.atm_types[cidx] = mds.Atom(sigma=cur_sig,
                                                epsilon=cur_eps,
                                                energy_unit="kcal/mol")

            # add old key as key and new key as value
            atm_key_old_new[cidx+1] = cidx

        # ityp:      atom type as number (= atom key)
        # imass:     mass of atom type
        # itypamb:   amber force field name of atom type
        # atypemass: dictionary with atom type masses and amber names assigned
        #            to atom-keys

        # type-id (key), mass, type-name
        atypemass = {}
        for ityp, imass, itypamb in zip(section_atom_type_index,
                                        section_mass,
                                        section_amber_atom_type):
            # iterate through names- and mass-list, overwrite duplicates
            # -> only single atom-types remain
            atypemass[ityp] = [imass, itypamb]

        for akey in atypemass:
            akey_new = atm_key_old_new[akey]  # translate old key
            self.atm_types[akey_new].weigh   = atypemass[akey][0]
            self.atm_types[akey_new].sitnam  = atypemass[akey][1]

        # /// bond-types
        bnd_key_old_new = {}
        for cbnd_id, (bnd_fconst, bnd_r0) in enumerate(zip(section_bond_force_constant,
                                                           section_bond_equil_value)):
            self.bnd_types[cbnd_id] = mds.Bond(prm1=bnd_fconst,
                                               prm2=bnd_r0,
                                               energy_unit="kcal/mol")
            bnd_key_old_new[cbnd_id+1] = cbnd_id

        # /// angle-types
        ang_key_old_new = {}
        for cang_id, (ang_fconst, ang_r0) in enumerate(zip(section_angle_force_constant,
                                                           section_angle_equil_value)):
            self.ang_types[cang_id] = mds.Angle(prm1=ang_fconst,
                                                prm2=ang_r0,
                                                energy_unit="kcal/mol",
                                                angle_unit="rad")
            ang_key_old_new[cang_id+1] = cang_id

        # /// dihedral- and improper-types
        #     assign dihedrals and impropers to "dih" and "imp"
        dih_imp_dict = agphf.unmask_imp(section_dihedrals_inc_hydrogen,
                                        section_dihedrals_without_hydrogen)

        # dicts where we can look up, which old dihedral-/improper-key points to the new key
        dih_old_new = {}
        imp_old_new = {}

        dih_cntr = 0  # dihedral key-counter, for consecutive numbering
        imp_cntr = 0  # improper key-counter, for consecutive numbering

        dih_key_old_new = {}
        #imp_key_old_new = {}
        for cdih_id, (dih_fconst, dih_prd, dih_phase) in enumerate(zip(section_dihedral_force_constant,
                                                                       section_dihedral_periodicity,
                                                                       section_dihedral_phase)):
            cur_key = cdih_id+1
            dih_prd = int(dih_prd)

            if dih_imp_dict[cur_key] == "dih":
                self.dih_types[dih_cntr] = mds.Dihedral(prm_k=dih_fconst,
                                                        prm_n=dih_prd,
                                                        prm_d=dih_phase,
                                                        energy_unit="kcal/mol",
                                                        angle_unit="rad")
                # pointer new-key -> old-key
                dih_old_new[cur_key] = dih_cntr
                dih_key_old_new[dih_cntr+1] = dih_cntr
                dih_cntr += 1  # right enumeration for dih-key (starting at 1)

            elif dih_imp_dict[cur_key] == "imp":
                self.imp_types[imp_cntr] = mds.Improper(prm_k=dih_fconst,
                                                        prm_n=dih_prd,
                                                        prm_d=dih_phase,
                                                        energy_unit="kcal/mol",
                                                        angle_unit="rad")
                # check if prm_k and prm_n fulfill cvff convention (see lammps)
                agphf.check_cvff_compatibility(dih_phase, dih_prd)
                # pointer new-key -> old-key
                imp_old_new[cur_key] = imp_cntr
                #imp_key_old_new[imp_cntr+1] = imp_cntr
                imp_cntr += 1  # right enumeration for imp-key (starting at 1)
            else:
                raise RuntimeError("***Undefined dihedral! Something went totally wrong!")

        # /// arrange data - topology section ///
        # /// atoms
        atm_id_old_new = {}
        for iatmkey, (itype, ichge, isitnam) in enumerate(zip(section_atom_type_index,
                                                              section_charge,
                                                              section_atom_name)):
            # since amber is/was based on fortran, increment atom-ids by 1
            self.atoms.append(mds.Atom(atm_id=iatmkey,
                                       atm_key=atm_key_old_new[itype],
                                       chge=ichge,
                                       sitnam=isitnam,
                                       grp_id=0))

            # save into the corresponding dictionary to translate amber-ids to internal ids
            #self.atm_idx_id[iatmkey] = iatmkey+1
            #self.atm_id_idx[iatmkey+1] = iatmkey
            atm_id_old_new[iatmkey+1] = iatmkey

        # gather residue information
        for cur_atm_nr in range(nres):
            cur_res = section_residue_label[cur_atm_nr]  # residue name

            if cur_atm_nr == 0:  # first
                pptr = 0
            else:
                pptr = section_residue_pointer[cur_atm_nr]  # present atom-pointer

            try:
                # atom-pointer after present atom-pointer
                iptr = section_residue_pointer[cur_atm_nr+1]
            except(IndexError):
                iptr = None  # None, if no pointer after present pointer

            for iatm in self.atoms[pptr:iptr]:
                iatm.res = cur_res

        # /// bonds
        self.bonds = agphf.relocate_parsed(section_bonds_inc_hydrogen,
                                           self.bonds,
                                           "bonds",
                                           bnd_key_old_new,
                                           atm_id_old_new)

        self.bonds = agphf.relocate_parsed(section_bonds_without_hydrogen,
                                           self.bonds,
                                           "bonds",
                                           bnd_key_old_new,
                                           atm_id_old_new)

        # /// angles
        self.angles = agphf.relocate_parsed(section_angles_inc_hydrogen,
                                            self.angles,
                                            "angles",
                                            ang_key_old_new,
                                            atm_id_old_new)

        self.angles = agphf.relocate_parsed(section_angles_without_hydrogen,
                                            self.angles,
                                            "angles",
                                            ang_key_old_new,
                                            atm_id_old_new)

        # /// dihedrals
        self.dihedrals = agphf.relocate_parsed(section_dihedrals_inc_hydrogen,
                                               self.dihedrals,
                                               "dihedrals",
                                               dih_old_new,
                                               atm_id_old_new)

        self.dihedrals = agphf.relocate_parsed(section_dihedrals_without_hydrogen,
                                               self.dihedrals,
                                               "dihedrals",
                                               dih_old_new,
                                               atm_id_old_new)

        # /// impropers
        self.impropers = agphf.relocate_parsed(section_dihedrals_inc_hydrogen,
                                               self.impropers,
                                               "impropers",
                                               imp_old_new,
                                               atm_id_old_new)

        self.impropers = agphf.relocate_parsed(section_dihedrals_without_hydrogen,
                                               self.impropers,
                                               "impropers",
                                               imp_old_new,
                                               atm_id_old_new)

        # delete variables not needed (at the moment)
        del prmtop_version
        del prmtop_title
        del numextra

        # if no ncopy is specified
        try:
            del ncopy
        except UnboundLocalError:
            pass

        # delete unneeded variables (at least at the moment there is no need)
        try:
            del section_number_excluded_atoms
            del section_scee_scale_factor
            del section_scnb_scale_factor
            del section_solty
            del section_radii
        except UnboundLocalError:
            pass

        # convert amber charges to normal ones
        print("***Prmtop-Info: Converting prmtop-charges.")
        self._convert_prmtop_charges()

        # convert impropers to match cvff (i.e. cos of prm_d)
        print("***Prmtop-Info: Converting impropers to cvff-style.")
        self._amb_imp2cvff()

        # form molecules by given bonds
        self.fetch_molecules_by_bonds()

    def _convert_prmtop_charges(self):
        """
        Charges are multiplied by 18.2223 (reasons unknown).
        """
        for cur_atm in self.atoms:
            cur_atm.chge /= 18.2223

    def _amb_imp2cvff(self):
        """
        Convert prm_d of improper to rad (to match cvff convention).
        """
        for iidx in self.imp_types:
            self.imp_types[iidx].cvff_prm_d()

    def read_inpcrd(self, inpcrd):
        """
        Read the contents of the amber inpcrd-file. Time and Temp will not
        be read!
        """
        print("***Inpcrd-Info: Reading Amber-Coordinates-File!")
        tmp_ts = []

        with open(inpcrd, "r") as inpcrd_in:
            title = inpcrd_in.readline()
            line = next(inpcrd_in)
            coords = int(line.split()[0])

            # number of following lines (6 coordinate entries per line)
            num_coords_lines = int(math.ceil(coords/2))

            # read coordinates line by line
            for _ in range(num_coords_lines):
                line = next(inpcrd_in)
                line_coords = [float(i) for i in line.split()]
                # define xyz-coordinates per line
                coords1 = np.array(line_coords[0:3])
                coords2 = np.array(line_coords[3:])
                tmp_ts.append(coords1)

                if coords2 != []:
                    tmp_ts.append(coords2)

            # check if coordinates match self.atoms
            if self.atoms != [] and len(self.atoms) != len(tmp_ts):
                raise Warning("Number of coordinates does not match number of atoms!")

            # add coordinates as first frame
            self.ts_coords.append(tmp_ts)

        del title

    def read_mdcrd(self, mdcrd):
        """
        Read an amber trajectory file (currently only coordinates can be read).
        """
        #TODO UNTESTED, PROBABLY WILL NOT WORK (RIGHT)
        print("***Mdcrd-Info: Reading Amber-Coordinates-File!")
        num_atms = len(self.atoms)
        num_coords = num_atms*3

        with open(mdcrd, "r") as mdcrd_in:
            reading = True
            title = next(mdcrd_in)

            # number of lines that describe one frame
            lines_per_frame = math.ceil(num_coords/10)
            cur_coords = []

            while reading is True:
                # read lines with coordinates
                for _ in range(lines_per_frame):
                    line = mdcrd_in.next().split()
                    coords = [float(i) for i in line]
                    cur_coords.append(coords)

                # coords per frame
                cur_frame  = cur_coords[:num_coords]

                # split array into sub-arrays
                cur_frame = np.array(cur_frame)
                cur_frame = np.split(cur_frame, 3)
                self.ts_coords.append(cur_frame)

                # remaining coords which are part of the next frame
                cur_coords = cur_coords[num_coords:]

        del title

    def write_mdcrd(self, mdcrd, *frame_ids):
        """
        Write an amber trajectory file.
        Source: http://ambermd.org/formats.html#trajectory
        """
        print("***Mdcrd-Info: Writing Amber-Coordinates-File!")
        flatten_coords = []

        with open(mdcrd, "w") as mdcrd_out:
            title = "Written output coordinates.\n"
            mdcrd_out.write(title)

            for frame_id in frame_ids:
                # flatten list
                flatten_coords = [coord for
                                  xyz_coords in self.ts_coords[frame_id] for
                                  coord in xyz_coords]

                for cidx, ccoord in enumerate(flatten_coords):
                    # write coordinates
                    mdcrd_out.write("{:>8.3f}".format(ccoord))

                    # newline after 10 coordinates
                    cidx += 1
                    if cidx % 10 == 0:
                        mdcrd_out.write("\n")
