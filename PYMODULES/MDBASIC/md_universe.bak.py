from __future__ import print_function, division
#import numpy as np
from natsort import natsorted
import copy
import md_universe_helper_functions as mduh
import md_linked_cells as mdlc

__version__ = "2017-06-21"


class Universe(object):
    """
    Generate lists for stuff almost all MD-Programs have in common, such as lists for:
        > atom-types
        > bond-types
        > angle-types
        > dihedral-types
    """
    def __init__(self):
        self.atoms       = []  # instances of Atom(); general stuff
        self.bonds       = []  # instances of Bond()
        self.angles      = []  # instances of Angle()
        self.dihedrals   = []  # instances of Dihedral()
        self.impropers   = []  # instances of Improper()
        self.atm_types   = {}  # instances of Atom(); force field stuff
        self.bnd_types   = {}  # instances of Bond(); force field stuff
        self.ang_types   = {}  # instances of Angle(); force field stuff
        self.dih_types   = {}  # instances of Dihedral(); force field stuff
        self.imp_types   = {}  # instances of Improper(); force field stuff
        self.long_range  = []  # instances of LongRange(); force field stuff
        self.ts_coords   = []  # all coordinates of all frames
        self.ts_boxes    = []  # instances of Box() of each frame
        self.ts_lnk_cls  = []  # instances of LinkedCells() of each frame
        #self.ts_steps   = []  # deprecate?
        self.molecules   = []  # containing lists with idxs of atoms that form a molecule
        self.atm_id_idx  = {}  # key: atom-id, value: atom-index (can be the same)
        self.atm_idx_id  = {}  # key: atom-index, value: atom-id (can be the same)

    def _atm_bnd_ang_dih_imp_id_idx(self, entry=None, vice_versa=False):
        """
        Convert all lammps-indices (ids) to internal indices (idx) and vice versa.

        Input:
            > entry             str; Container which atom-ids shall be converted:
                                atms|bnds|angs|dihs|imps
            > vice_versa        boolean; convert internal indices to given (lammps-)indices
        """
        if vice_versa:
            _dict = self.atm_idx_id
        else:
            _dict = self.atm_id_idx

        if entry == "atms":
            for cur_atm_idx, cur_atm in enumerate(self.atoms):
                    self.atoms[cur_atm_idx].atm_id = _dict[cur_atm.atm_id]
        elif entry == "bnds":
            for cur_bnd_idx, cur_bnd in enumerate(self.bonds):
                self.bonds[cur_bnd_idx].atm_id1 = _dict[cur_bnd.atm_id1]
                self.bonds[cur_bnd_idx].atm_id2 = _dict[cur_bnd.atm_id2]
        elif entry == "angs":
            for cur_ang_idx, cur_ang in enumerate(self.angles):
                self.angles[cur_ang_idx].atm_id1 = _dict[cur_ang.atm_id1]
                self.angles[cur_ang_idx].atm_id2 = _dict[cur_ang.atm_id2]
                self.angles[cur_ang_idx].atm_id3 = _dict[cur_ang.atm_id3]
        elif entry == "dihs":
            for cur_dih_idx, cur_dih in enumerate(self.dihedrals):
                self.dihedrals[cur_dih_idx].atm_id1 = _dict[cur_dih.atm_id1]
                self.dihedrals[cur_dih_idx].atm_id2 = _dict[cur_dih.atm_id2]
                self.dihedrals[cur_dih_idx].atm_id3 = _dict[cur_dih.atm_id3]
                self.dihedrals[cur_dih_idx].atm_id4 = _dict[cur_dih.atm_id4]
        elif entry == "imps":
            for cur_imp_idx, cur_imp in enumerate(self.impropers):
                self.impropers[cur_imp_idx].atm_id1 = _dict[cur_imp.atm_id1]
                self.impropers[cur_imp_idx].atm_id2 = _dict[cur_imp.atm_id2]
                self.impropers[cur_imp_idx].atm_id3 = _dict[cur_imp.atm_id3]
                self.impropers[cur_imp_idx].atm_id4 = _dict[cur_imp.atm_id4]
        else:
            raise ValueError("Container of name {} unknown!".format(entry))

    def create_linked_cells(self, frame_id=0, rcut_a=2, rcut_b=2, rcut_c=2):
        """
        Categorize atoms by their membership in certain linked cells.
        Create those cells and affiliate the atoms there.

        Input:
            > frame_id      int;
            > ra            float; factor to divide side length a by.
            > rb            float;              -"-             b by.
            > rc            float;              -"-             c by.
        """
        # get lattice box vectors (if cell is not already fractional)
        this_boxtype = self.ts_boxes[frame_id].boxtype
        tmp_copy_box = copy.copy(self.ts_boxes[frame_id])

        # convert (the copy) of the current box-type to a fractional box-type
        this_coords_type = None
        if this_boxtype == "lammps":
            tmp_copy_box.box_lmp2lat()
            this_coords_type = "cartesian"
        elif this_boxtype == "cartesian":
            tmp_copy_box.box_cart2lat()
            this_coords_type = "cartesian"
        else:
            this_coords_type = "lattice"

        # do some init magic (convert cartesian coordinates, etc.)
        linked_cells = mdlc.LinkedCells(self.ts_coords[frame_id],
                                        tmp_copy_box.ltc_a,
                                        tmp_copy_box.ltc_b,
                                        tmp_copy_box.ltc_c,
                                        tmp_copy_box.ltc_alpha,
                                        tmp_copy_box.ltc_beta,
                                        tmp_copy_box.ltc_gamma,
                                        coords_type=this_coords_type)

        linked_cells.create_lnk_cells(rcut_a, rcut_b, rcut_c)
        self.ts_linked_cells.append(linked_cells)

    def _merge_dicts(self, dict2, atm_types=True, bnd_types=False, ang_types=False,
                     dih_types=False, imp_types=False):
        """
        Merge two entries (self and another) of atm_types, bnd_types,
        ang_types, dih_types or imp_types.
        The already existing dict (e.g. self.atm_types) will be readjusted during
        the process!

        Input:
            > dict2     dict; 2nd atm_types, 2nd bnd_types, etc. to load on the
                        already existing ones. If types are not already present,
                        append them to the existing ones, otherwise ignore them.
            > atm_types, bnd_types, etc.    boolean; tells which kind of type
                                            is being processed.
        """
        print("HI")

        if atm_types:
            ref_dict = mduh.merge_entries(self.atm_types, dict2)
        elif bnd_types:
            ref_dict = mduh.merge_entries(self.bnd_types, dict2)
        elif ang_types:
            ref_dict = mduh.merge_entries(self.ang_types, dict2)
        elif dih_types:
            ref_dict = mduh.merge_entries(self.dih_types, dict2)
        elif imp_types:
            ref_dict = mduh.merge_entries(self.imp_types, dict2)
        else:
            raise Warning("Nothing to merge!")

        print(ref_dict)
        return ref_dict

    def _add_atoms(self, atoms2, atm_types2=True):
        """
        Add instances of class Atom to already existing list self.atoms
        """

        # merge both atom-type-dictionaries and append those types which have
        # no overlap with the existing ones
        if atm_types2 is True:
            ref_atm_type = self._merge_dicts(atm_types2, atm_types=True)

            # change keys so they fit the new atom-type-references
            for catm in atoms2:
                catm.atm_key = ref_atm_type[catm.atm_key]

        # append new atoms (with right keys)
        self.atoms += atoms2

    def _add_bonds(self, bnd_types2, bonds2):
        ref_bnd_type = self._merge_dicts(bonds2, atm_types=True)
        print(ref_bnd_type)

    def _sort(self, keyword, atm_keyword=None):
        """
        Sort
            > atoms by id, key, charge, x-, y- or z-coordinates
            > bonds by atom-id 1
            > angles by atom-id 1
            > dihedrals by atom-id 1
            > impropers by atom-id 1
        """
        # sort atoms by id or atm_keyword
        if keyword == "atoms":

            # standard setting = sorting by atom-id
            if atm_keyword == "atm_key":
                self.atoms = natsorted(self.atoms, key=mduh.get_atm_key)
            elif atm_keyword == "atm_grp":
                self.atoms = natsorted(self.atoms, key=mduh.get_atm_atm_grp)
            elif atm_keyword == "sitnam":
                self.atoms = natsorted(self.atoms, key=mduh.get_atm_sitnam)
            elif atm_keyword == "res":
                self.atoms = natsorted(self.atoms, key=mduh.get_atm_res)
            elif atm_keyword == "chge":
                self.atoms = natsorted(self.atoms, key=mduh.get_atm_chge)
            else:
                self.atoms = natsorted(self.atoms, key=mduh.get_atm_id)

        # sort bonds, angles, dihedrals, impropers by first atom-id
        elif keyword == "bonds":
            self.bonds = natsorted(self.bonds, key=mduh.get_atm_id1)
        elif keyword == "angles":
            self.angles = natsorted(self.angles, key=mduh.get_atm_id1)
        elif keyword == "dihedrals":
            self.dihedrals = natsorted(self.dihedrals, key=mduh.get_atm_id1)
        elif keyword == "impropers":
            self.impropers = natsorted(self.impropers, key=mduh.get_atm_id1)
        else:
            pass

    def refresh(self):
        """
        Renumbering/reordering of atoms (i.e. sort), e.g. after molecules
        or atoms were deleted during a simulation.
        """
        # dict where we can look up previous atom-indices of new atom-indices
        assigned_atm_ids = {}

        # /// refresh atoms  ///
        self._sort("atoms")  # presort atoms by id
        for catm_id, iatm in enumerate(self.atoms):
            # reassign atom-ids == overwrite old ones
            assigned_atm_ids[iatm.atm_id] = catm_id + 1
            iatm.atm_id = catm_id + 1

        # /// refresh bonds ///
        for cbnd_id, ibnd in enumerate(self.bonds):
            # reassign atom-ids
            ibnd.atm_id1 = assigned_atm_ids[ibnd.atm_id1]
            ibnd.atm_id2 = assigned_atm_ids[ibnd.atm_id2]
        self._sort("bonds")

        # /// refresh angles ///
        for cang_id, iang in enumerate(self.angles):
            # reassign atom-ids
            iang.atm_id1 = assigned_atm_ids[iang.atm_id1]
            iang.atm_id2 = assigned_atm_ids[iang.atm_id2]
            iang.atm_id3 = assigned_atm_ids[iang.atm_id3]
            iang.ang_id  = cang_id + 1
        self._sort("angles")

        # /// refresh dihedrals ///
        for cdih_id, idih in enumerate(self.dihedrals):
            # reassign atom-ids
            idih.atm_id1 = assigned_atm_ids[idih.atm_id1]
            idih.atm_id2 = assigned_atm_ids[idih.atm_id2]
            idih.atm_id3 = assigned_atm_ids[idih.atm_id3]
            idih.atm_id4 = assigned_atm_ids[idih.atm_id4]
            idih.dih_id  = cdih_id + 1
        self._sort("dihedrals")

        # /// refresh impropers
        for cimp_id, iimp in enumerate(self.impropers):
            # reassign atom-ids
            iimp.atm_id1 = assigned_atm_ids[iimp.atm_id1]
            iimp.atm_id2 = assigned_atm_ids[iimp.atm_id2]
            iimp.atm_id3 = assigned_atm_ids[iimp.atm_id3]
            iimp.atm_id4 = assigned_atm_ids[iimp.atm_id4]
            iimp.imp_id = cimp_id + 1
        self._sort("impropers")

        # refresh indices
        self._atoms_idx()

    def _atoms_idx(self):
        """
        Get indices for each atom-id. Helper function for refresh-method.
        Return:
            atm_id_idx      dict; keys = atom-ids, values = atom-indices
        """
        for cidx, catm in enumerate(self.atoms):
            self.atm_id_idx[catm.atm_id] = cidx

        # create dict which is vice versa to atm_id_idx
        self.atm_idx_id = {y: x for x, y in self.atm_id_idx.iteritems()}

    def delete_atoms(self, *atm_idxs):
        """
        Delete a single atom (and all its angles, bonds, dihedrals, impropers)
        optional: delete whole molecule, the atom belongs to.
        atm_idxs = atom-ids to delete (getting more during search for molecules)
        molecule:   'yes' or something else
        """
        # fetch whole molecules
        atoms2delete_idxs = atm_idxs

        # /// delete further entries with chosen atom-ids
        for catm_idx in atoms2delete_idxs:
            self.atoms[catm_idx] = None

        # convert idxs from atoms_to_delete to ids since the following entries
        # only contain atom-ids (not idxs); note: sets are fast for searching!
        atoms2delete_ids = set([self.atm_idx_id[i] for i in atoms2delete_idxs])

        # crawl bonds
        for bnd_idx, cbnd in enumerate(self.bonds):
            if cbnd.atm_id1 in atoms2delete_ids or cbnd.atm_id2 in atoms2delete_ids:
                self.bonds[bnd_idx] = None

        # crawl angles
        for ang_idx, cang in enumerate(self.angles):
            if (cang.atm_id1 in atoms2delete_ids or cang.atm_id2 in atoms2delete_ids or
               cang.atm_id3 in atoms2delete_ids):
                self.angles[ang_idx] = None

        # crawl dihedrals
        for dih_idx, cdih in enumerate(self.dihedrals):
            if (cdih.atm_id1 in atoms2delete_ids or cdih.atm_id2 in atoms2delete_ids or
               cdih.atm_id3 in atoms2delete_ids or cdih.atm_id4 in atoms2delete_ids):
                self.dihedrals[dih_idx] = None

        # crawl impropers
        for imp_idx, cimp in enumerate(self.impropers):
            if (cimp.atm_id1 in atoms2delete_ids or cimp.atm_id2 in atoms2delete_ids or
               cimp.atm_id3 in atoms2delete_ids or cimp.atm_id4 in atoms2delete_ids):
                self.impropers[imp_idx] = None

        # delete Nones from lists (quite fast)
        self.atoms = mduh._del_nones(self.atoms)
        self.bonds = mduh._del_nones(self.bonds)
        self.angles = mduh._del_nones(self.angles)
        self.dihedrals = mduh._del_nones(self.dihedrals)
        self.impropers = mduh._del_nones(self.impropers)
        self.molecules = mduh._del_nones(self.molecules)

    def ui_convert_units(self, energy_unit_out='eV', ang_unit_out=False,
                         cvff_style=False):
        """
        Convert read units of force constants, angles, charges to given format.
        Standard AMBER-Units are kcal/mol*A**-2 (force constants)
        Conversion is directed from Attribute 'energy_unit' of each entry!
        energy unit from which is converted is already defined by each
        mds-class (Atom, Bond, Angle, etc.) during instantiating
        (should be at least).
        energy_unit_out:    str; 'eV'|'kcal/mol'|'kj/mol'
        ang_unit_out:       str; 'deg'|'rad'
        Sources:    http://lammps.sandia.gov/doc/units.html
        """
        if self.bnd_types:
            for cur_bndtype in self.bnd_types:
                self.bnd_types[cur_bndtype].convert_energy_unit(energy_unit_out)

        if self.ang_types:
            for cur_angtype in self.ang_types:
                self.ang_types[cur_angtype].convert_energy_unit(energy_unit_out)
                if ang_unit_out:
                    self.ang_types[cur_angtype].convert_angle_unit(ang_unit_out)

        if self.dih_types:
            for cur_dihtype in self.dih_types:
                self.dih_types[cur_dihtype].convert_energy_unit(energy_unit_out)
                if ang_unit_out:
                    self.dih_types[cur_dihtype].convert_angle_unit(ang_unit_out)

        if self.imp_types:
            for cur_imptype in self.imp_types:
                self.imp_types[cur_imptype].convert_energy_unit(energy_unit_out)
                if cvff_style is True:
                    self.imp_types[cur_imptype].cvff_prm_d()

    def combine_universes(self, universe2, frame_id, mode="append"):
        """
        Combine the data if two instances (self, universe2) in different modes
        such as append, complement or overwrite.
        """
        #TODO CONCATENATE BOXES AND ATOMIC COORDINATES
        #TODO CONCATENATE MOLECULES, REASSIGN ATOM-IDS
        # appending mode (do not replace existing entries, append everything
        # from universe 2)
        if mode == "append":

            # atoms-, bonds-, angles-, dihedrals-, impropers-section
            for entry in ("atoms", "bonds", "angles", "dihedrals", "impropers"):

                # define u1_ce (current entry) and u1_cet (current entry type - force field info)
                if entry == "atoms":
                    u1_ce = self.atoms
                    u1_cet = self.atm_types
                    u2_ce = universe2.atoms
                    u2_cet = universe2.atm_types
                    # translation between former and new atom-ids
                    u1_atm_id_old_new = {}
                    u2_atm_id_old_new = {}
                elif entry == "bonds":
                    u1_ce = self.bonds
                    u1_cet = self.bnd_types
                    u2_ce = universe2.bonds
                    u2_cet = universe2.bnd_types
                elif entry == "angles":
                    u1_ce = self.angles
                    u1_cet = self.ang_types
                    u2_ce = universe2.angles
                    u2_cet = universe2.ang_types
                elif entry == "dihedrals":
                    u1_ce = self.dihedrals
                    u1_cet = self.dih_types
                    u2_ce = universe2.dihedrals
                    u2_cet = universe2.dih_types
                else:
                    u1_ce = self.impropers
                    u1_cet = self.imp_types
                    u2_ce = universe2.impropers
                    u2_cet = universe2.imp_types

                # disassemble former keys and values (force field) - universe 1
                u1_len_ce = len(u1_ce)
                u1_cet_tp_keys = u1_cet.keys()
                u1_cet_tp_vals = u1_cet.values()
                u1_len_cet_tp_keys = len(u1_cet_tp_keys)

                # disassemble former keys and values (force field) - universe 2
                u2_len_ce = len(u2_ce)
                u2_cet_tp_keys = u2_cet.keys()
                u2_cet_tp_vals = u2_cet.values()

                # for translation between former and new keys (both universes)
                u1_cet_id_old_new = {}
                u2_cet_id_old_new = {}

                # redefine indices
                all_cet_tp_keys = range(u1_len_cet_tp_keys+len(u2_cet_tp_keys))

                # concatenate values of current entry
                all_cet_tp_vals = u1_cet_tp_vals + u2_cet_tp_vals

                # translate between old and new entry-keys (universe 1)
                u1_cet_keys_old_new = zip(u1_cet_tp_keys, all_cet_tp_keys[:u1_len_cet_tp_keys])
                u1_cet_keys_old_new = dict(u1_cet_keys_old_new)

                # translate between old and new entry-keys (universe 2)
                u2_cet_keys_old_new = zip(u2_cet_tp_keys, all_cet_tp_keys[u1_len_cet_tp_keys:])
                u2_cet_keys_old_new = dict(u2_cet_keys_old_new)

                # concatenation to new entry for all keys of current entry
                u1_cet = dict(zip(all_cet_tp_keys, all_cet_tp_vals))

                # bonds, angles, dihedrals, impropers - concatenation
                if entry == "atoms":
                    u1_ce.extend(universe2.atoms)
                elif entry == "bonds":
                    u1_ce.extend(universe2.bonds)
                elif entry == "angles":
                    u1_ce.extend(universe2.angles)
                elif entry == "dihedrals":
                    u1_ce.extend(universe2.dihedrals)
                else:
                    u1_ce.extend(universe2.impropers)

                # renew atoms-, bonds-, angles-, dihedrals- and impropers-entries
                for cimp_id in xrange(u1_len_ce+u2_len_ce):

                    # translate force-field-keys and atom-ids of universe 1 part
                    if cimp_id < u1_len_ce:

                        # translate old force field keys to new keys
                        if entry == "atoms":
                            u1_old_cet_key = u1_ce[cimp_id].atm_key
                        elif entry == "bonds":
                            u1_old_cet_key = u1_ce[cimp_id].bnd_key
                        elif entry == "angles":
                            u1_old_cet_key = u1_ce[cimp_id].ang_key
                        elif entry == "dihedrals":
                            u1_old_cet_key = u1_ce[cimp_id].dih_key
                        else:
                            u1_old_cet_key = u1_ce[cimp_id].imp_key

                        u1_new_cet_key = u1_cet_keys_old_new[u1_old_cet_key]

                        if entry == "atoms":
                            u1_ce[cimp_id].atm_key = u1_new_cet_key
                            old_ce_id = u1_ce[cimp_id].atm_id
                        elif entry == "bonds":
                            u1_ce[cimp_id].bnd_key = u1_new_cet_key
                            old_ce_id = u1_ce[cimp_id].bnd_id
                        elif entry == "angles":
                            u1_ce[cimp_id].ang_key = u1_new_cet_key
                            old_ce_id = u1_ce[cimp_id].ang_id
                        elif entry == "dihedrals":
                            u1_ce[cimp_id].dih_key = u1_new_cet_key
                            old_ce_id = u1_ce[cimp_id].dih_id
                        else:
                            u1_ce[cimp_id].imp_key = u1_new_cet_key
                            old_ce_id = u1_ce[cimp_id].imp_id

                        # differ definition between atom-ids and other ids since we need
                        # to translate between old and new atom-ids for bonds,
                        # angles, dihedrals and impropers but not vice versa
                        if entry == "atoms":
                            u1_atm_id_old_new[old_ce_id] = cimp_id
                        else:
                            u1_cet_id_old_new[old_ce_id] = cimp_id

                        # change atom-ids for bonds, angles, dihedrals and
                        # impropers - universe 1

                        if entry != "atoms":
                            # atom-id-1
                            old_atm_id1 = u1_ce[cimp_id].atm_id1
                            new_atm_id1 = u1_atm_id_old_new[old_atm_id1]
                            u1_ce[cimp_id].atm_id1 = new_atm_id1

                            # atom-id-2
                            old_atm_id2 = u1_ce[cimp_id].atm_id2
                            new_atm_id2 = u1_atm_id_old_new[old_atm_id2]
                            u1_ce[cimp_id].atm_id2 = new_atm_id2

                        # atom-id-3
                        if entry in ("angles", "dihedrals", "impropers"):
                            old_atm_id3 = u1_ce[cimp_id].atm_id3
                            new_atm_id3 = u1_atm_id_old_new[old_atm_id3]
                            u1_ce[cimp_id].atm_id3 = new_atm_id3

                        # atom-id-4
                        if entry in ("dihedrals", "impropers"):
                            old_atm_id4 = u1_ce[cimp_id].atm_id4
                            new_atm_id4 = u1_atm_id_old_new[old_atm_id4]
                            u1_ce[cimp_id].atm_id4 = new_atm_id4

                    # translate force-field-keys and atom-ids of universe 2 part
                    else:

                        # get old force-field-key (universe 2)
                        if entry == "atoms":
                            u2_old_cet_key = u1_ce[cimp_id].atm_key
                        elif entry == "bonds":
                            u2_old_cet_key = u1_ce[cimp_id].bnd_key
                        elif entry == "angles":
                            u2_old_cet_key = u1_ce[cimp_id].ang_key
                        elif entry == "dihedrals":
                            u2_old_cet_key = u1_ce[cimp_id].dih_key
                        else:
                            u2_old_cet_key = u1_ce[cimp_id].imp_key

                        u2_new_cet_key = u2_cet_keys_old_new[u2_old_cet_key]

                        # redefine force-field-key (universe 2)
                        if entry == "atoms":
                            u1_ce[cimp_id].atm_key = u2_new_cet_key
                            old_ce_id = u1_ce[cimp_id].atm_id
                        elif entry == "bonds":
                            u1_ce[cimp_id].bnd_key = u2_new_cet_key
                            old_ce_id = u1_ce[cimp_id].bnd_id
                        elif entry == "angles":
                            u1_ce[cimp_id].ang_key = u2_new_cet_key
                            old_ce_id = u1_ce[cimp_id].ang_id
                        elif entry == "dihedrals":
                            u1_ce[cimp_id].dih_key = u2_new_cet_key
                            old_ce_id = u1_ce[cimp_id].dih_id
                        else:
                            u1_ce[cimp_id].imp_key = u2_new_cet_key
                            old_ce_id = u1_ce[cimp_id].imp_id

                        # change atom-ids for bonds, angles, dihedrals and
                        # impropers - universe 2
                        if entry == "atoms":
                            # assign new atom-ids to each atom
                            u2_atm_id_old_new[old_ce_id] = cimp_id
                        else:
                            u2_cet_id_old_new[old_ce_id] = cimp_id

                        if entry != "atoms":
                            # atom-id-1
                            old_atm_id1 = u1_ce[cimp_id].atm_id1
                            new_atm_id1 = u2_atm_id_old_new[old_atm_id1]
                            u1_ce[cimp_id].atm_id1 = new_atm_id1
                            # atom-id-2
                            old_atm_id2 = u1_ce[cimp_id].atm_id2
                            new_atm_id2 = u2_atm_id_old_new[old_atm_id2]
                            u1_ce[cimp_id].atm_id2 = new_atm_id2

                        # atom-id-3
                        if entry in ("angles", "dihedrals", "impropers"):
                            old_atm_id3 = u1_ce[cimp_id].atm_id3
                            new_atm_id3 = u2_atm_id_old_new[old_atm_id3]
                            u1_ce[cimp_id].atm_id3 = new_atm_id3

                        # atom-id-4
                        if entry in ("dihedrals", "impropers"):
                            old_atm_id4 = u1_ce[cimp_id].atm_id4
                            new_atm_id4 = u2_atm_id_old_new[old_atm_id4]
                            u1_ce[cimp_id].atm_id4 = new_atm_id4

                    # assign new ids
                    if entry == "atoms":
                        u1_ce[cimp_id].atm_id = cimp_id
                        self.atm_types = u1_cet
                    elif entry == "bonds":
                        u1_ce[cimp_id].bnd_id = cimp_id
                        self.bnd_types = u1_cet
                    elif entry == "angles":
                        u1_ce[cimp_id].ang_id = cimp_id
                        self.ang_types = u1_cet
                    elif entry == "dihedrals":
                        u1_ce[cimp_id].dih_id = cimp_id
                        self.dih_types = u1_cet
                    else:
                        u1_ce[cimp_id].imp_id = cimp_id
                        self.imp_types = u1_cet

        # merge universe 1 and universe 2 (where possible)
        elif mode == "merge":
            # atoms-, bonds-, angles-, dihedrals-, impropers-section
            for entry in ("atoms", "bonds", "angles", "dihedrals", "impropers"):

                # define u1_ce (current entry) and u1_cet (current entry type - force field info)
                if entry == "atoms":
                    u1_ce = self.atoms
                    u1_cet = self.atm_types
                    u2_ce = universe2.atoms
                    u2_cet = universe2.atm_types
                elif entry == "bonds":
                    u1_ce = self.bonds
                    u1_cet = self.bnd_types
                    u2_ce = universe2.bonds
                    u2_cet = universe2.bnd_types
                elif entry == "angles":
                    u1_ce = self.angles
                    u1_cet = self.ang_types
                    u2_ce = universe2.angles
                    u2_cet = universe2.ang_types
                elif entry == "dihedrals":
                    u1_ce = self.dihedrals
                    u1_cet = self.dih_types
                    u2_ce = universe2.dihedrals
                    u2_cet = universe2.dih_types
                else:
                    u1_ce = self.impropers
                    u1_cet = self.imp_types
                    u2_ce = universe2.impropers
                    u2_cet = universe2.imp_types

                # merge entries, get dict for u2 to translate between old and new keys
                u2_cet_old_new_tmp = mduh.merge_entries(u1_cet, u2_cet)

                # renumber keys and merge new keys with left over values from
                # u1 and u2 (all_new_entry)
                u1_cet_keys = u1_cet.keys()  # old keys
                all_new_keys = range(len(u1_cet_keys))  # new keys
                all_new_vals = u1_cet.values()
                all_new_entry = dict(zip(all_new_keys, all_new_vals))

                # get dictionaries to translate between original keys from
                # u1, u2 and new keys from all_new_entry
                u1_cet_old_new = {}
                u2_cet_old_new = {}

                # create an intermediate dictionary to translate between
                # newest and renewed keys
                d_newer_newest = dict(zip(u1_cet_keys, all_new_keys))

                for key in u1_cet.keys():
                    u1_cet_old_new[key] = d_newer_newest[key]

                for old_key, new_key in zip(u2_cet_old_new_tmp.keys(),
                                            u2_cet_old_new_tmp.values()):
                    u2_cet_old_new[old_key] = d_newer_newest[new_key]

                # renew ff-section
                if entry == "atoms":
                    self.atm_types = all_new_entry
                elif entry == "bonds":
                    self.bnd_types = all_new_entry
                elif entry == "angles":
                    self.ang_types = all_new_entry
                elif entry == "dihedrals":
                    self.dih_types = all_new_entry
                else:
                    self.imp_types = all_new_entry

                # renew other section
                #TODO RENEW FF-KEYS AND ATM-IDS
                for i in xrange(len(u1_ce)+len(u2_ce)):
                    pass
                break
