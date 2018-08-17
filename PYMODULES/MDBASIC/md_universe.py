"""
Module which combines classes from submodules.

MD-Universe is a class that has everything a system needs to be useful for
carrying out molecular dynamics. Given bond information it groups atoms to
molecules.
"""

from __future__ import print_function, division
import numpy as np
import copy
import itertools as it
import time
from natsort import natsorted
import Transformations as cgt
import ag_geometry as agm
import ag_cryst as agc
import ag_vectalg as agv
import md_stars as mds
import md_elements as mde
import md_box as mdb
import md_linked_cells as mdlc
import md_universe_helper_functions as mduh
import networkx
from networkx.algorithms.components.connected import connected_components
import rmsd
import pdb

__version__ = "2018-01-10"

#=================
# Helper Functions
#=================


def to_graph(l):
    G = networkx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(to_edges(part))
    return G


def to_edges(l):
    """
    treat `l` as a Graph and returns it's edges
    to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current


class Universe(object):
    """
    Generate lists for stuff almost all MD-Programs have in common, such as lists for:
        > atom-types
        > bond-types
        > angle-types
        > dihedral-types
    """
    def __init__(self):
        """
        Initialize all containers which will later be needed (maybe).
        """
        # atom-property section
        self.atm_types   = {}  # instances of Atom(); force field stuff
        self.atoms       = []  # instances of Atom(); general stuff
        # molecular topology sections
        self.bonds       = []  # instances of Bond()
        self.angles      = []  # instances of Angle()
        self.dihedrals   = []  # instances of Dihedral()
        self.impropers   = []  # instances of Improper()
        self.molecules   = []  # containing lists with idxs of atoms that form a molecule
        # force field sections
        self.bnd_types   = {}  # instances of Bond(); force field stuff
        self.ang_types   = {}  # instances of Angle(); force field stuff
        self.dih_types   = {}  # instances of Dihedral(); force field stuff
        self.imp_types   = {}  # instances of Improper(); force field stuff
        self.pair_types  = []  # holds all pair-coefficients
        # coordinate and box sections
        self.ts_coords   = []  # all coordinates of all frames
        self.ts_forces   = []  # all forces of all frames
        #self.ts_velocs   = []  # all velocities of all frames
        self.ts_boxes    = []  # instances of Box() of each frame
        self.ts_lnk_cls  = []  # instances of LinkedCells() of each frame

    # COMMON-STUFF -------------------------------------------------------------------
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
            assigned_atm_ids[iatm.atm_id] = catm_id
            iatm.atm_id = catm_id

        # /// refresh bonds ///
        self._sort("bonds")
        for cbnd_id, ibnd in enumerate(self.bonds):
            # reassign atom-ids
            ibnd.atm_id1 = assigned_atm_ids[ibnd.atm_id1]
            ibnd.atm_id2 = assigned_atm_ids[ibnd.atm_id2]
            ibnd.bnd_id = cbnd_id

        # /// refresh angles ///
        self._sort("angles")
        for cang_id, iang in enumerate(self.angles):
            # reassign atom-ids
            iang.atm_id1 = assigned_atm_ids[iang.atm_id1]
            iang.atm_id2 = assigned_atm_ids[iang.atm_id2]
            iang.atm_id3 = assigned_atm_ids[iang.atm_id3]
            iang.ang_id  = cang_id

        # /// refresh dihedrals ///
        self._sort("dihedrals")
        for cdih_id, idih in enumerate(self.dihedrals):
            # reassign atom-ids
            idih.atm_id1 = assigned_atm_ids[idih.atm_id1]
            idih.atm_id2 = assigned_atm_ids[idih.atm_id2]
            idih.atm_id3 = assigned_atm_ids[idih.atm_id3]
            idih.atm_id4 = assigned_atm_ids[idih.atm_id4]
            idih.dih_id  = cdih_id

        # /// refresh impropers
        self._sort("impropers")
        for cimp_id, iimp in enumerate(self.impropers):
            # reassign atom-ids
            iimp.atm_id1 = assigned_atm_ids[iimp.atm_id1]
            iimp.atm_id2 = assigned_atm_ids[iimp.atm_id2]
            iimp.atm_id3 = assigned_atm_ids[iimp.atm_id3]
            iimp.atm_id4 = assigned_atm_ids[iimp.atm_id4]
            iimp.imp_id = cimp_id

        # /// refresh molecules

    def delete_atoms(self, *atoms2delete):
        """
        Delete a single atom (and all its angles, bonds, dihedrals, impropers)
        optional: delete whole molecule, the atom belongs to.
        atm_idxs = atom-ids to delete (getting more during search for molecules)
        molecule:   'yes' or something else

        Input:
            > same_molecule     boolean; choose if the whole molecule shall be
                                deleted; needs to be enabled except no bonds
                                were defined!
            > atoms2delete      list; atom indices that shall be deleted
        """
        # delete whole molecules
        atoms2delete = set(self.same_molecule_as(False, *atoms2delete))

        # /// delete entries
        for frame_id in xrange(len(self.ts_coords)):
            self.ts_coords[frame_id] = [coord for coord_idx, coord in enumerate(self.ts_coords[frame_id]) if
                                        coord_idx not in atoms2delete]

        self.atoms     = [atm for idx, atm in enumerate(self.atoms) if
                          idx not in atoms2delete]
        self.bonds     = [bnd for bnd in self.bonds if
                          bnd.atm_id1 not in atoms2delete and
                          bnd.atm_id2 not in atoms2delete]
        self.angles    = [ang for ang in self.angles if
                          ang.atm_id1 not in atoms2delete and
                          ang.atm_id2 not in atoms2delete and
                          ang.atm_id3 not in atoms2delete]
        self.dihedrals = [dih for dih in self.dihedrals if
                          dih.atm_id1 not in atoms2delete and
                          dih.atm_id2 not in atoms2delete and
                          dih.atm_id3 not in atoms2delete and
                          dih.atm_id4 not in atoms2delete]
        self.impropers = [imp for imp in self.impropers if
                          imp.atm_id1 not in atoms2delete and
                          imp.atm_id2 not in atoms2delete and
                          imp.atm_id3 not in atoms2delete and
                          imp.atm_id4 not in atoms2delete]

        # delete molecules
        for molecule_idx, molecule in enumerate(self.molecules):
            for atom in molecule:
                if atom in atoms2delete:
                    self.molecules[molecule_idx] = None
                    break

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
        if self.atm_types:
            for cur_atmtype in self.atm_types:
                self.atm_types[cur_atmtype].convert_energy_unit(energy_unit_out)

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

        if self.pair_types:
            print("***UI-Convert-Warning: Conversion for pair_types not implemented yet (but for atoms it is)!")
            pass

    def mix_pair_types(self, mode="ii", mix_style="arithmetic", to_file=None, debug=False):
        """
        Calculate sigma_ij and epsilon_ij by utilizing the Lorentz-Berthelot rules
        Sources:    https://en.wikipedia.org/wiki/Combining_rules
                    DOI: 10.1002/andp.18812480110
                    https://www.researchgate.net/post/Where_are_the_charge_sigma_epsilon_parameters_in_GAFF
                    http://lammps.sandia.gov/doc/pair_modify.html
        vdw1, vdw2 are instances of class VdW.

        mix arithmetic:
            sigma_ij = (sigma_i + sigma_j) / 2
            epsilon_ij = sqrt(epsilon_i * epsilon_j)
        mix geometric:
            sigma_ij = sqrt(sigma_i * sigma_j)
            epsilon_ij = sqrt(epsilon_i * epsilon_j)
        to_file:
            write mixing to a file with filename being the value of 'to_file'
        """
        # delete previous values
        if debug is True:
            print("***Mixing Pair Types Info: Deleting previous pair_types!")

        if mode == "ii":
            if self.pair_types != []:
                self.pair_types = []

            for i in self.atm_types:
                # define current atom-type
                atm_i = self.atm_types[i]

                if hasattr(atm_i, "sigma")and hasattr(atm_i, "epsilon"):
                    # mix vdw of current atom-type with itself
                    sigma_ii, epsilon_ii = atm_i.mix_ij(atm_i.sigma,
                                                        atm_i.epsilon,
                                                        mix=mix_style)
                    # create instance of LongRange()
                    cii = mds.LongRange(lr_key="lj",
                                        atm_key_i=i,
                                        atm_key_j=i,
                                        sigma_ij=sigma_ii,
                                        epsilon_ij=epsilon_ii,
                                        pairs="ii"
                                        )

                    # write results to file
                    if to_file is not None:
                        with open(to_file, "a") as pair_file:
                            pair_file.write("{:<5}{:<5}{:>10}{:>10}\n".format(cij.atm_key_i, cij.atm_key_j, cij.sigma_ij, cij.epsilon_ij))
                    else:
                        self.pair_types.append(cii)

        elif mode == "ij":
            #TODO Mix and do not overwrite existing styles, i.e. i and j must be
            #TODO renamed when appending takes place between different files

            for i, j in it.combinations_with_replacement(self.atm_types, 2):
                atm_i = self.atm_types[i]
                atm_j = self.atm_types[j]

                if hasattr(atm_i, "sigma")and hasattr(atm_i, "epsilon") \
                   and hasattr(atm_j, "sigma")and hasattr(atm_j, "epsilon"):
                    sigma_ij, epsilon_ij = atm_i.mix_ij(atm_j.sigma,
                                                        atm_j.epsilon,
                                                        mix=mix_style)
                    # create instance of LongRange()
                    cij = mds.LongRange(lr_key="lj",
                                        atm_key_i=i,
                                        atm_key_j=j,
                                        sigma_ij=sigma_ij,
                                        epsilon_ij=epsilon_ij,
                                        pairs="ij")

                    # write results to file
                    if to_file is not None:
                        with open(to_file, "a") as pair_file:
                            pair_file.write("{:<5}{:<5}{:<20}{:<20}\n".format(cij.atm_key_i, cij.atm_key_j, cij.epsilon_ij, cij.sigma_ij))
                    else:
                        self.pair_types.append(cij)

        else:
            raise IOError("Wrong mode. Allowed: ii|ij")

    def extend_universe(self,
                        universe2,
                        u1_frame_id=0,
                        u2_frame_id=0,
                        mode="merge"):
        #TODO adapt box from 1 or 2 argument
        #TODO ATTN THIS STUFF HERE NEEDS REVISION TO STAY IN LINE WITH THE OTHER
        #TODO METHODS. EVERYTHING SHOULD BE DONE IN A MUCH CLEARER WAY!
        """
        Combine the data if two instances (self, universe2_copy) in different modes
        such as append, complement or overwrite.
        Input:
            > universe2     instance of mdu.Universe; second universe to append,
                            merge or complement with this instance of Universe
            > u1_frame_id   int; frame-id of frame from universe 1 to merge with u2
            > u2_frame_id   int; frame-id of frame from universe 2 to merge with u1
            > mode          str; 'merge', 'append' or 'complement'
                            The way how the given instance shall be extended
        """
        #TODO write pair coeffs with right parameters (ii is enough, in lammps fix mix paircoeffs arithmetic)
        # appending mode (do not replace existing entries, append everything
        # from universe 2)
        # make a copy since we do not want to change the original universe_2
        universe2_copy = copy.deepcopy(universe2)

        # atoms-, bonds-, angles-, dihedrals-, impropers-section
        for entry in ("atoms", "bonds", "angles", "dihedrals", "impropers"):

            # define u1_ce (current entry) and u1_ffe (current entry type - force field info)
            if entry == "atoms":
                u1_ce = self.atoms
                u1_ffe = self.atm_types
                u2_ce = universe2_copy.atoms
                u2_ffe = universe2_copy.atm_types
                # translation between former and new atom-ids
                u1_atm_id_old_new = {}
                u2_atm_id_old_new = {}
            elif entry == "bonds":
                u1_ce = self.bonds
                u1_ffe = self.bnd_types
                u2_ce = universe2_copy.bonds
                u2_ffe = universe2_copy.bnd_types
            elif entry == "angles":
                u1_ce = self.angles
                u1_ffe = self.ang_types
                u2_ce = universe2_copy.angles
                u2_ffe = universe2_copy.ang_types
            elif entry == "dihedrals":
                u1_ce = self.dihedrals
                u1_ffe = self.dih_types
                u2_ce = universe2_copy.dihedrals
                u2_ffe = universe2_copy.dih_types
            else:
                u1_ce = self.impropers
                u1_ffe = self.imp_types
                u2_ce = universe2_copy.impropers
                u2_ffe = universe2_copy.imp_types

            # disassemble former keys and values (force field) - universe 1
            u1_len_ce = len(u1_ce)
            u1_ffe_keys = u1_ffe.keys()
            u1_ffe_vals = u1_ffe.values()
            u1_len_ffe_keys = len(u1_ffe_keys)

            # disassemble former keys and values (force field) - universe 2
            u2_len_ce = len(u2_ce)
            u2_ffe_keys = u2_ffe.keys()
            u2_ffe_vals = u2_ffe.values()
            u1_len_ffe_keys = len(u2_ffe_vals)

            # for translation between former and new keys (both universes)
            u1_ffe_keys_old_new = {}
            u2_ffe_keys_old_new = {}

            if mode in ("merge", "append"):
                if mode == "append":
                    # redefine indices
                    combined_ffe_keys = range(u1_len_ffe_keys+len(u2_ffe_keys))

                    # concatenate values of current entry
                    combined_ffe_vals = u1_ffe_vals + u2_ffe_vals

                    # translate between old and new entry-keys (universe 1)
                    u1_ffe_keys_old_new = zip(u1_ffe_keys, combined_ffe_keys[:u1_len_ffe_keys])
                    u1_ffe_keys_old_new = dict(u1_ffe_keys_old_new)

                    # translate between old and new entry-keys (universe 2)
                    u2_ffe_keys_old_new = zip(u2_ffe_keys, combined_ffe_keys[u1_len_ffe_keys:])
                    u2_ffe_keys_old_new = dict(u2_ffe_keys_old_new)

                # merge universe 1 and universe 2 (where possible)
                elif mode == "merge":
                    # check similarities in current force-field-entry and create
                    # a dictionary to translate between old ff-keys and new ones
                    combined_ffe_keys = u1_ffe_keys[:]
                    combined_ffe_vals = u1_ffe_vals[:]
                    u1_ffe_keys_old_new = dict(zip(u1_ffe_keys, u1_ffe_keys))
                    # attributes of class instances of force-field-entry from u2
                    # as dictionary
                    for u2_cffe_key, u2_cvalue in zip(u2_ffe_keys, u2_ffe_vals):
                        # dictionary with attributes and their values for u1
                        u2_cattr = dict(u2_cvalue.__iter__())
                        u2_cattr_keys = u2_cattr.keys()
                        u2_cattr_vals = u2_cattr.values()

                        # attributes of class instances of force-field-entry from u1
                        # as dictionary
                        for u1_cffe_key, u1_cvalue in zip(u1_ffe_keys, u1_ffe_vals):
                            same = None
                            # dictionary with attributes and their values for u1
                            u1_cattr = dict(u1_cvalue.__iter__())
                            u1_cattr_keys = u1_cattr.keys()
                            u1_cattr_vals = u1_cattr.values()

                            # check if attributes are the same
                            if u2_cattr_keys != u1_cattr_keys:
                                same = False

                            # check if attributes' values are the same
                            if u2_cattr_vals == u1_cattr_vals:
                                same = True
                                break
                            else:
                                same  = False

                        if same is False:
                            combined_ffe_keys.append(len(combined_ffe_keys)+1)
                            combined_ffe_vals.append(u2_cvalue)
                            u2_ffe_keys_old_new[u2_cffe_key] = len(combined_ffe_keys)
                        else:
                            u2_ffe_keys_old_new[u2_cffe_key] = u1_cffe_key

                    # make keys start with key-id 0
                    # create dict to translate between old and new keys
                    combined_ffe_keys_tmp = {}
                    for idx, ckey in enumerate(combined_ffe_keys):
                        combined_ffe_keys_tmp[ckey] = idx

                    # dict to translate between old u1-ff-keys and new ones
                    for ckey in u1_ffe_keys_old_new:
                        u1_ffe_keys_old_new[ckey] = combined_ffe_keys_tmp[ckey]

                    # dict to translate between old u2-ff-keys and new ones
                    for ckey in u2_ffe_keys_old_new:
                        u2_ffe_keys_old_new[ckey] = combined_ffe_keys_tmp[u2_ffe_keys_old_new[ckey]]

                    # redefine keys for force field entry to start with 0
                    combined_ffe_keys = range(len(combined_ffe_keys))

                else:
                    raise TypeError("Wrong keyword for argument 'mode'")

                # concatenation to new entry for all keys of current entry
                u1_ffe = dict(zip(combined_ffe_keys, combined_ffe_vals))

                # bonds, angles, dihedrals, impropers - concatenation
                u1_ce.extend(u2_ce)

                # renew atoms-, bonds-, angles-, dihedrals- and impropers-entries
                for ce_idx in xrange(len(u1_ce)):

                    # translate force-field-keys and atom-ids of universe 1 part
                    if ce_idx < u1_len_ce:

                        # translate old force field keys to new keys
                        if entry == "atoms":
                            u1_old_ffe_key = u1_ce[ce_idx].atm_key
                        elif entry == "bonds":
                            u1_old_ffe_key = u1_ce[ce_idx].bnd_key
                        elif entry == "angles":
                            u1_old_ffe_key = u1_ce[ce_idx].ang_key
                        elif entry == "dihedrals":
                            u1_old_ffe_key = u1_ce[ce_idx].dih_key
                        else:
                            u1_old_ffe_key = u1_ce[ce_idx].imp_key

                        # translate force field key of current atom, bond, etc.
                        # to newly defined one
                        u1_new_ffe_key = u1_ffe_keys_old_new[u1_old_ffe_key]

                        if entry == "atoms":
                            u1_ce[ce_idx].atm_key = u1_new_ffe_key
                            old_ce_id = u1_ce[ce_idx].atm_id
                        elif entry == "bonds":
                            u1_ce[ce_idx].bnd_key = u1_new_ffe_key
                            old_ce_id = u1_ce[ce_idx].bnd_id
                        elif entry == "angles":
                            u1_ce[ce_idx].ang_key = u1_new_ffe_key
                            old_ce_id = u1_ce[ce_idx].ang_id
                        elif entry == "dihedrals":
                            u1_ce[ce_idx].dih_key = u1_new_ffe_key
                            old_ce_id = u1_ce[ce_idx].dih_id
                        else:
                            u1_ce[ce_idx].imp_key = u1_new_ffe_key
                            old_ce_id = u1_ce[ce_idx].imp_id

                        # differ definition between atom-ids and other ids since we need
                        # to translate between old and new atom-ids for bonds,
                        # angles, dihedrals and impropers but not vice versa
                        if entry == "atoms":
                            u1_atm_id_old_new[old_ce_id] = ce_idx

                        # change atom-ids for bonds, angles, dihedrals and
                        # impropers - universe 1

                        if entry != "atoms":
                            # atom-id-1
                            old_atm_id1 = u1_ce[ce_idx].atm_id1
                            new_atm_id1 = u1_atm_id_old_new[old_atm_id1]
                            u1_ce[ce_idx].atm_id1 = new_atm_id1

                            # atom-id-2
                            old_atm_id2 = u1_ce[ce_idx].atm_id2
                            new_atm_id2 = u1_atm_id_old_new[old_atm_id2]
                            u1_ce[ce_idx].atm_id2 = new_atm_id2

                        # atom-id-3
                        if entry in ("angles", "dihedrals", "impropers"):
                            old_atm_id3 = u1_ce[ce_idx].atm_id3
                            new_atm_id3 = u1_atm_id_old_new[old_atm_id3]
                            u1_ce[ce_idx].atm_id3 = new_atm_id3

                        # atom-id-4
                        if entry in ("dihedrals", "impropers"):
                            old_atm_id4 = u1_ce[ce_idx].atm_id4
                            new_atm_id4 = u1_atm_id_old_new[old_atm_id4]
                            u1_ce[ce_idx].atm_id4 = new_atm_id4

                    # translate force-field-keys and atom-ids of universe 2 part
                    else:

                        # get old force-field-key (universe 2)
                        if entry == "atoms":
                            u2_old_ffe_key = u1_ce[ce_idx].atm_key
                        elif entry == "bonds":
                            u2_old_ffe_key = u1_ce[ce_idx].bnd_key
                        elif entry == "angles":
                            u2_old_ffe_key = u1_ce[ce_idx].ang_key
                        elif entry == "dihedrals":
                            u2_old_ffe_key = u1_ce[ce_idx].dih_key
                        else:
                            u2_old_ffe_key = u1_ce[ce_idx].imp_key

                        u2_new_ffe_key = u2_ffe_keys_old_new[u2_old_ffe_key]

                        # redefine force-field-key (universe 2)
                        if entry == "atoms":
                            u1_ce[ce_idx].atm_key = u2_new_ffe_key
                            old_ce_id = u1_ce[ce_idx].atm_id
                        elif entry == "bonds":
                            u1_ce[ce_idx].bnd_key = u2_new_ffe_key
                            old_ce_id = u1_ce[ce_idx].bnd_id
                        elif entry == "angles":
                            u1_ce[ce_idx].ang_key = u2_new_ffe_key
                            old_ce_id = u1_ce[ce_idx].ang_id
                        elif entry == "dihedrals":
                            u1_ce[ce_idx].dih_key = u2_new_ffe_key
                            old_ce_id = u1_ce[ce_idx].dih_id
                        else:
                            u1_ce[ce_idx].imp_key = u2_new_ffe_key
                            old_ce_id = u1_ce[ce_idx].imp_id

                        # change atom-ids for bonds, angles, dihedrals and
                        # impropers - universe 2
                        if entry == "atoms":
                            # assign new atom-ids to each atom
                            u2_atm_id_old_new[old_ce_id] = ce_idx

                        if entry != "atoms":
                            # atom-id-1
                            old_atm_id1 = u1_ce[ce_idx].atm_id1
                            new_atm_id1 = u2_atm_id_old_new[old_atm_id1]
                            u1_ce[ce_idx].atm_id1 = new_atm_id1
                            # atom-id-2
                            old_atm_id2 = u1_ce[ce_idx].atm_id2
                            new_atm_id2 = u2_atm_id_old_new[old_atm_id2]
                            u1_ce[ce_idx].atm_id2 = new_atm_id2

                        # atom-id-3
                        if entry in ("angles", "dihedrals", "impropers"):
                            old_atm_id3 = u1_ce[ce_idx].atm_id3
                            new_atm_id3 = u2_atm_id_old_new[old_atm_id3]
                            u1_ce[ce_idx].atm_id3 = new_atm_id3

                        # atom-id-4
                        if entry in ("dihedrals", "impropers"):
                            old_atm_id4 = u1_ce[ce_idx].atm_id4
                            new_atm_id4 = u2_atm_id_old_new[old_atm_id4]
                            u1_ce[ce_idx].atm_id4 = new_atm_id4

                    # assign new ids
                    if entry == "atoms":
                        u1_ce[ce_idx].atm_id = ce_idx
                        self.atm_types = u1_ffe
                    elif entry == "bonds":
                        u1_ce[ce_idx].bnd_id = ce_idx
                        self.bnd_types = u1_ffe
                    elif entry == "angles":
                        u1_ce[ce_idx].ang_id = ce_idx
                        self.ang_types = u1_ffe
                    elif entry == "dihedrals":
                        u1_ce[ce_idx].dih_id = ce_idx
                        self.dih_types = u1_ffe
                    else:
                        u1_ce[ce_idx].imp_id = ce_idx
                        self.imp_types = u1_ffe

            else:  # mode is complement; mostly for merging two lammps files; does not work yet
                # let us check if the user knows what he is doing
                if u1_len_ce != u2_len_ce and (u1_len_ce != 0 and u2_len_ce != 0):
                    raise Warning("Entries of both universes have different sizes! No way to complement this")

                if u1_len_ffe_keys != u1_len_ffe_keys and (u1_len_ffe_keys != 0 and u1_len_ffe_keys != 0):
                    raise Warning("Force Field entries of both universes have different sizes! No way to complement this")

                    # attributes of class instances of force-field-entry from u2
                    # as dictionary
                    for u1_cffe_key, u1_cvalue, u2_cffe_key, u2_cvalue in zip(
                            zip(u2_ffe_keys, u2_ffe_vals), zip(u1_ffe_keys, u1_ffe_vals)):
                        # dictionary with attributes and their values
                        u1_cattr = dict(u1_cvalue.__iter__())
                        u1_cattr_keys = u1_cattr.keys()
                        u1_cattr_vals = u1_cattr.values()
                        u2_cattr = dict(u2_cvalue.__iter__())
                        u2_cattr_keys = u2_cattr.keys()
                        u2_cattr_vals = u2_cattr.values()

        if mode in ("merge", "append"):
            # do not do anything with boxes, it is up to the user to adjust the box
            # size
            print("***Warning: Ignoring Box-Info from universe 2!")

            # append newly formed molecules
            u1_len_molecules = len(self.molecules)
            self.molecules.extend(universe2_copy.molecules)

            for cmol_idx in xrange(len(self.molecules)):

                # adjust atom-ids of molecules from universe 1
                if cmol_idx < u1_len_molecules:
                    self.molecules[cmol_idx] = [u1_atm_id_old_new[i] for i in
                                                self.molecules[cmol_idx]]
                # adjust atom-ids of molecules from universe 2
                else:
                    self.molecules[cmol_idx] = [u2_atm_id_old_new[i] for i in
                                                self.molecules[cmol_idx]]

            # append coordinates
            if universe2_copy.ts_coords != []:

                if self.ts_coords == []:
                    self.ts_coords = np.array(self.ts_coords)

                self.ts_coords[u1_frame_id] = np.concatenate(
                    (self.ts_coords[u1_frame_id], universe2_copy.ts_coords[u2_frame_id]),
                    axis=0)

        #TODO merge stuff here needs urgent revision since it was coded quickly and dirty
        #TODO this here is just a temporary fix for pair types stuff
        self.pair_types = []  # erase all old pair coefficients
        self.mix_pair_types(mode="ii", mix_style="arithmetic", debug=False)

    # COORDINATE-STUFF ------------------------------------------------------------------------------
    def _append_4th_dimension(self, frame_id):
        """
        Helper method to append a 4th dimension to a whole frame.
        """
        num_atms = len(self.ts_coords[frame_id])
        filling_ones = np.ones((num_atms, 1))
        self.ts_coords[frame_id] = np.concatenate((self.ts_coords[frame_id],
                                                   filling_ones), axis=1)

    def mm_atm_coords(self, frame_id, Mx, copy, *atom_ids):
        """
        Returns a copy matrix-multiplied coordinates, i.e. does not alter
        the actual coordinates of Universe.

        frame_id    int; id of frame to manipulate the coordinates
        Mx          np-array; matrix to multiply the coordinates with
        copy        boolean; return the altered coordinates as a copy (True) or
                    change the internal coordinates (False)
        atom_ids    int or list/tuple; ids of atoms to alter
        """
        copy_coords = np.copy(self.ts_coords[frame_id])

        # Check matrix' columns
        if len(Mx[0]) == 4:
            num_atms = len(copy_coords)
            filling_ones = np.ones((num_atms, 1))
            copy_coords = np.concatenate((copy_coords, filling_ones), axis=1)

        for cidx in atom_ids:
            ccoords = np.matmul(Mx, copy_coords[cidx])
            copy_coords[cidx] = ccoords

        # delete 4th dimension after coordinates were altered
        if len(Mx[0]) == 4:
            copy_coords = np.delete(copy_coords, 3, axis=1)

        if copy is True:
            return copy_coords
        else:
            self.ts_coords[frame_id] = copy_coords

    def get_com(self, frame_id, *atom_ids):
        """
        Calculate the center of mass for given atoms (by their atom ids).
        frame       int or None; index of frame
        atom_idx    boolean; interpret atom_ids as indices (True) or as
                    atom-ids (False)
        atm_ids     ints or list; all atom indices or -ids that should be
                    considered for the com calculation
                    if using a list for atom_ids, do no forget the asterisk!
                    e.g. *mylist
        """
        masses = []
        coords = []
        # chosen frame
        cur_ts_coords = self.ts_coords[frame_id]

        # translate atom-id to atom-index (if atm_idx is False)
        for cidx in atom_ids:
            # get mass of current atom
            ckey  = self.atoms[cidx].atm_key
            cmass = self.atm_types[ckey].weigh
            masses.append(cmass)

            # get coordinates of current atom
            ccoords = cur_ts_coords[cidx]
            coords.append(ccoords)

        # calculate com
        center = agm.get_com(coords, masses)
        return center

    def get_cog(self, frame_id, *atom_ids):
        """
        Calculate the center of geometry for given atoms (by their atom ids).
        frame       int or None; index of frame
        atom_idx    boolean; interpret atom_ids as indices (True) or as
                    atom-ids (False)
        atm_ids     ints or list; all atom indices or -ids that should be
                    considered for the com calculation
                    if using a list for atom_ids, do no forget the asterisk!
                    e.g. *mylist
        """
        coords = []
        # chosen frame
        cur_ts_coords = self.ts_coords[frame_id]

        # translate atom-id to atom-index (if atm_idx is False)
        for cidx in atom_ids:
            # get coordinates of current atom
            ccoords = cur_ts_coords[cidx]
            coords.append(ccoords)

        # calculate com
        center = agm.get_cog(coords)
        return center

    def def_boxes_by_coords(self):
        """
        If no box vectors are given, define (a) rectangular box(es) by
        given atomic coordinates.
        """

        for cid, cframe in enumerate(self.ts_coords):
            # (re)define initial coordinates
            x_smallest, y_smallest, z_smallest = 1e20, 1e20, 1e20
            x_biggest, y_biggest, z_biggest = -1e20, -1e20, -1e20

            for ccoord in cframe:

                if ccoord[0] < x_smallest:
                    x_smallest = ccoord[0]

                if ccoord[1] < y_smallest:
                    y_smallest = ccoord[1]

                if ccoord[2] < z_smallest:
                    z_smallest = ccoord[2]

                if ccoord[0] > x_biggest:
                    x_biggest = ccoord[0]

                if ccoord[1] > y_biggest:
                    y_biggest = ccoord[1]

                if ccoord[2] > z_biggest:
                    z_biggest = ccoord[2]

            ca = abs(x_smallest-x_biggest)
            cb = abs(y_smallest-y_biggest)
            cc = abs(z_smallest-z_biggest)
            cbox = mdb.Box(boxtype="lattice",
                           ltc_a=ca,
                           ltc_b=cb,
                           ltc_c=cc,
                           ltc_alpha=np.radians(90),
                           ltc_beta=np.radians(90),
                           ltc_gamma=np.radians(90))

            # try to modify a current given box or append if none is given
            try:
                self.ts_boxes[cid] = cbox
            except IndexError:
                self.ts_boxes.append(cbox)

    def get_system_radius(self, frame_id):
        """
        Get the distance between the center of geometry (of the molecule) to that
        atom that is furthest away from it.

        Input:
            > frame_id      int;

        Returns:
            > radius        float; biggest distance from center of geometry
            > cog           np-array; center of geometry

        """
        # center of geometry
        cog = self.get_cog(0, *range(len(self.atoms)))
        radius = 0  # very small number for initial distance

        # get distances for all atoms between center of geometry and current atom
        for ccoord in self.ts_coords[frame_id]:
            cur_vt = cog - ccoord
            cur_dist = np.linalg.norm(cur_vt)

            if cur_dist > radius:
                radius = cur_dist

        return radius

    def get_mol_radius(self, frame_id, *atm_ids):
        """
        Get radius of a molecule.

        Input:
            > frame_id      int;

        Returns:
            > radius        float; biggest distance from center of geometry
            > cog           np-array; center of geometry
        """
        cog = self.get_cog(frame_id, *atm_ids)
        radius = 0

        for cidx in atm_ids:
            cur_vt = self.ts_coords[frame_id][cidx] - cog
            cur_dist = np.linalg.norm(cur_vt)

            if cur_dist > radius:
                radius = cur_dist

        return (radius, cog)

    def create_linked_cells(self, frame_id, rcut_a=2, rcut_b=2, rcut_c=2):
        """
        Categorize atoms by their membership in certain linked cells.
        Create those cells and affiliate the atoms there.

        Input:
            > frame_id      int;
            > ra            float; side length of sub-cell vector a
            > rb            float; side length of sub-cell vector b
            > rc            float; side length of sub-cell vector c
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
        self.ts_lnk_cls.append(linked_cells)

    def chk_atm_dist(self,
                     frame_id=-1,
                     min_dist=0.80,
                     exclude_same_molecule=True,
                     get_aggregates=False,
                     debug=False):
        """
        Check the inter atomic distances inside the sub cells and between
        atoms of neighboring sub cells. Each sub cell has 26 neighbor cells.

        Input:
            > min_dist      float; minimal distance between two atoms
            > periodic      boolean; check also atoms which lie in periodic boxes.
            > exclude_same_molecule boolean; do not compare atoms of the same molecule
            > get_aggregates boolean; check if all atoms are part of the same aggregate,
                                return an array of atom-indices if True for each aggregate

        Return:
            > close_contacts    list of lists; contains all atom-idx with closer
                                contacts as min_dist
        """
        if debug is True:
            # not sure if that is the case - i think this works with wrapped cells
            print("***Info: Atom distance check only works with unwrapped cells! " +
                  "Unwrap cells if necessary.")
        # PBC - get cartesian box vectors
        this_boxtype = self.ts_boxes[frame_id].boxtype
        tmp_box = copy.copy(self.ts_boxes[frame_id])

        # convert (the copy) of the current box-type to a cartesian box-type
        if this_boxtype == "lammps":
            tmp_box.box_lmp2cart()
        elif this_boxtype == "lattice":
            tmp_box.box_lat2cart()
        else:
            pass

        if exclude_same_molecule is True and debug is True:
            print("***Info: 'exclude_same_molecule' chosen. " +
                  "Groups should be assigned by molecule affiliation " +
                  "or this will fail!")

        # get maximum possible number of aggregates
        if get_aggregates is True:
            connected_groups = []

        # linked cells must have been defined before
        linked_cells = self.ts_lnk_cls[frame_id]

        # to check intra-sub-cellular distances convert the fractional coordinates
        # back to cartesian coordinates
        close_contacts = []

        # number n-1 of boxes to check (without the last one since we start counting from 0 to n-1)
        ra = linked_cells.ra
        rb = linked_cells.rb
        rc = linked_cells.rc

        if debug is True:
            print("***Info: Checking distances")
            start = time.time()

        # matrix cartesian to fractional (inverse of fractional matrix)
        M_fc = agc.M_fract2cart(linked_cells.ltc_a,
                                linked_cells.ltc_b,
                                linked_cells.ltc_c,
                                linked_cells.ltc_alpha,
                                linked_cells.ltc_beta,
                                linked_cells.ltc_gamma)

        for cra in xrange(ra):
            # check from a along box vector b
            for crb in xrange(rb):
                # finally check along box vector c (from a and b)
                for crc in xrange(rc):

                    if debug is True:
                        print("Checking current cell: ", cra, crb, crc)

                    # position of current cell
                    current_cell = linked_cells.linked_cells[cra][crb][crc]

                    # 0 = same cell, 1 = next cell along a
                    for i in [-1, 0, 1]:
                        cra_idx = cra + i

                        # 0 = same cell, 1 = next cell along b
                        for j in [-1, 0, 1]:
                            crb_idx = crb + j

                            # 0 = same cell, 1 = next cell along c
                            for k in [-1, 0, 1]:
                                crc_idx = crc + k

                                try:
                                    # check coordinates of current cell vs coordinates of cur_neighbor_cell
                                    cur_neighbor_cell = linked_cells.linked_cells[cra_idx][crb_idx][crc_idx]
                                except IndexError:
                                    continue

                                if debug is True:
                                    print("Checking neighbor cell: ", cra_idx, crb_idx, crc_idx)

                                for cidx_a in current_cell:
                                    for cidx_b in cur_neighbor_cell:
                                        # make a temporary copy of the current
                                        # coordinates since we do not want to
                                        # alter them permanently
                                        tmp_coords_cidx_b = np.copy(linked_cells.atm_coords[cidx_b])

                                        # skip atoms which are part of the same molecule (same grp_id)
                                        if exclude_same_molecule is True and (self.atoms[cidx_a].grp_id == self.atoms[cidx_b].grp_id):
                                            continue

                                        # shift atomic coordinates by one box-vector if beyond pbc
                                        if cra_idx == -1:
                                            tmp_coords_cidx_b[0] = linked_cells.atm_coords[cidx_b][0] - 1
                                        if crb_idx == -1:
                                            tmp_coords_cidx_b[1] = linked_cells.atm_coords[cidx_b][1] - 1

                                        if crc_idx == -1:
                                            tmp_coords_cidx_b[2] = linked_cells.atm_coords[cidx_b][2] - 1

                                        # get distance vector
                                        vt_a = linked_cells.atm_coords[cidx_a]
                                        vt_b = tmp_coords_cidx_b
                                        vt_ab = vt_b - vt_a
                                        # convert vector to cartesian coordinates
                                        vt_ab = np.matmul(M_fc, vt_ab)
                                        # get magnitude of distance vector
                                        cdist = np.linalg.norm(vt_ab)

                                        # compare if its magnitude is smaller than the given threshold
                                        if cdist != 0 and cdist <= min_dist:
                                            close_contacts.append(cidx_a)
                                            close_contacts.append(cidx_b)

                                            # add group id, but only if this combination is not already present
                                            if get_aggregates is True:
                                                if [self.atoms[cidx_a].grp_id, self.atoms[cidx_b].grp_id] not in connected_groups:
                                                    connected_groups.append([self.atoms[cidx_a].grp_id, self.atoms[cidx_b].grp_id])

                                        elif cdist == 0 and (cidx_a != cidx_b):
                                            if debug is True:
                                                print("***Warning: Distance between {} and {} is 0!".format(cidx_a, cidx_b))
                                            close_contacts.append(cidx_a)
                                            close_contacts.append(cidx_b)
                                        else:
                                            pass

        if debug is True:
            end = time.time()
            print("***Info: Distance search finished after: {} seconds.".format(end-start))

        # remove duplicates
        close_contacts = set(close_contacts)

        #pdb.set_trace()
        #if get_aggregates is True:
        #    print("***Info: Connected groups: ", connected_groups)

        # ==============================#
        # merge molecules to aggregates
        # ==============================#
        if get_aggregates is True:
            # create a dict with all groups
            connections = to_graph(connected_groups)
            aggregates = connected_components(connections)
            aggregates = [i for i in aggregates]  # generator to list
            #pdb.set_trace()
            return (close_contacts, aggregates)

        return close_contacts

    def unwrap_cell(self, frame_id=0):
        """
        Unwrap the unit cell, i.e. move all atoms to the quadrant that is positive
        in x, y and z direction.
        Caveat: Currently this only works with orthogonal unit cells.
        """
        print("***Info: Unwrapping cell")
        # get box vectors
        # PBC - get cartesian box vectors
        this_boxtype = self.ts_boxes[frame_id].boxtype
        cp_box = copy.copy(self.ts_boxes[frame_id])

        # convert box to fractional
        if this_boxtype == "cartesian":
            cp_box.box_cart2lat()
        elif this_boxtype == "lammps":
            cp_box.box_lmp2lat()
        else:
            pass

        # convert coordinates: cartesian -> fractional
        M_cf = agc.M_cart2fract(cp_box.ltc_a, cp_box.ltc_b, cp_box.ltc_c,
                                cp_box.ltc_alpha, cp_box.ltc_beta,
                                cp_box.ltc_gamma)

        # transpose coords for 3x3-matrix multiplication
        ts_coords_T = np.transpose(self.ts_coords[frame_id])
        matmul_coords_T = np.matmul(M_cf, ts_coords_T)
        # transpose back
        self.ts_coords[frame_id] = np.transpose(matmul_coords_T)

        # move all negative coordinates to the positive quadrant
        for idx in xrange(len(self.ts_coords[frame_id])):

            # transpose about a
            if self.ts_coords[frame_id][idx][0] < 0:
                self.ts_coords[frame_id][idx][0] += 1

            # transpose about b
            if self.ts_coords[frame_id][idx][1] < 0:
                self.ts_coords[frame_id][idx][1] += 1

            # transpose about c
            if self.ts_coords[frame_id][idx][2] < 0:
                self.ts_coords[frame_id][idx][2] += 1

        # convert coordinates: fractional -> cartesian
        M_cf = agc.M_fract2cart(cp_box.ltc_a, cp_box.ltc_b, cp_box.ltc_c,
                                cp_box.ltc_alpha, cp_box.ltc_beta,
                                cp_box.ltc_gamma)

        # transpose coords for 3x3-matrix multiplication
        ts_coords_T = np.transpose(self.ts_coords[frame_id])
        matmul_coords_T = np.matmul(M_cf, ts_coords_T)
        # transpose back
        self.ts_coords[frame_id] = np.transpose(matmul_coords_T)

    def replicate_cell(self, n_start=0, n_stop=0, direction="a", frame_id=-1, adjust_box=True):
        """
        Replicate the given cell n-times in direction of the crystallographic
        box-vectors a, b or c.
        """
        if n_stop - n_start != 0:
            print("***Info: Replicating cell")
            # replicate topology
            self.add_topology_replicate(abs(n_stop-n_start))

            # get box vectors
            # PBC - get cartesian box vectors
            this_boxtype = self.ts_boxes[frame_id].boxtype
            cp_box = copy.copy(self.ts_boxes[frame_id])

            # convert (the copy) of the current box-type to a fractional box-type
            if this_boxtype == "lammps":
                cp_box.box_lmp2cart()
            elif this_boxtype == "lattice":
                cp_box.box_lat2cart()
            else:
                pass

            # prepare coordinates for matrix magic
            cp_coords = np.copy(self.ts_coords[frame_id])
            num_atms = len(cp_coords)
            filling_ones = np.ones((num_atms, 1))
            cp_coords = np.concatenate((cp_coords, filling_ones), axis=1)
            del (num_atms, filling_ones)

            # transpose coords for matrix-multiplication
            cp_coords = np.transpose(cp_coords)
            shift_vt = None

            # only translation in positive direction -> starting from next unit
            # cell after current one
            if n_start == 0:
                scope = range(n_start+1, n_stop+1)
            else:
                scope = range(n_start, n_stop+1)

            for replica_id in scope:

                # skip id 0, as it is the original unit cell
                if replica_id == 0:
                    continue

                # calculate the vector to shift the coordinates
                if direction == "a":
                    shift_vt = np.array(cp_box.crt_a)
                elif direction == "b":
                    shift_vt = np.array(cp_box.crt_b)
                else:
                    shift_vt = np.array(cp_box.crt_c)

                # calculate total shift of the coordinates from the original unit cell
                total_shift_vt = shift_vt * replica_id

                # do some matrix-multiplication magic
                Tm = cgt.translation_matrix(total_shift_vt)
                cur_shifted_coords = np.matmul(Tm, cp_coords)
                cur_shifted_coords = np.transpose(cur_shifted_coords)
                cur_shifted_coords = np.delete(cur_shifted_coords, 3, axis=1)
                # add coordinates to existing ones
                self.ts_coords[frame_id] = np.vstack([self.ts_coords[frame_id], cur_shifted_coords])

            # adjust the box-size according to the multiplication of the vectors
            if adjust_box is True:
                # calculate the vector to shift the coordinates
                if direction == "a":
                    cp_box.crt_a = np.array(cp_box.crt_a) * abs(n_stop - n_start) + np.array(cp_box.crt_a)
                elif direction == "b":
                    cp_box.crt_b = np.array(cp_box.crt_b) * abs(n_stop - n_start) + np.array(cp_box.crt_b)
                else:
                    cp_box.crt_c = np.array(cp_box.crt_c) * abs(n_stop - n_start) + np.array(cp_box.crt_c)

                # convert (the copy) of the current box-type to a fractional box-type
                if this_boxtype == "lammps":
                    cp_box.box_cart2lmp()
                elif this_boxtype == "lattice":
                    cp_box.box_cart2lat()
                else:
                    pass

                self.ts_boxes[frame_id] = cp_box

    # MOLECULE-STUFF --------------------------------------------------------------------------
    def guess_atomtypes(self, by_mass=False, by_typename=False, overwrite=False):
        """
        Guess atom-names (sitnam) by mass. Overwrite existing ones!
        by_mass     boolean; guess atom name by given mass
        by_typename boolean; take atom name by given atom name from atom type
        overwrite   boolean; overwrite given atom names
        REMEMBER: THIS WORKS! IF NO SITNAM IS FOUND, CHECK IF CGCMM WAS PARSED
                  IN LAMMPS FILES!
        """
        # guess name by given mass
        if by_mass:
            for cur_atom in self.atoms:

                # overwrite entries having atom-names already assigned
                if overwrite:
                    ckey = cur_atom.atm_key
                    # get mass from atm_types-dictionary
                    mass = round(self.atm_types[ckey].weigh, 1)
                    cur_atom.sitnam = mde.elements[mass]

                # skip entries having atom-names already assigned
                else:

                    try:
                        if cur_atom.sitnam:
                            pass
                    # all entries without a sitnam
                    except AttributeError:

                        for atm_typ in self.atm_types:
                            ckey = cur_atom.atm_key
                            # get mass from atm_types-dictionary
                            mass = round(self.atm_types[ckey].weigh, 1)
                            cur_atom.sitnam = mde.elements[mass]

        # guess name by cgcmm-info from 'Masses'-section
        elif by_typename:
            for cur_atom in self.atoms:

                if overwrite:

                    try:
                        ckey = cur_atom.atm_key
                        cur_atom.sitnam = self.atm_types[ckey].sitnam
                    # no sitnam given
                    except AttributeError:
                        pass

                # do not overwrite existing entries
                else:

                    try:
                        # skip if sitnam is given
                        if cur_atom.sitnam:
                            pass
                    except AttributeError:
                        ckey = cur_atom.atm_key
                        try:
                            cur_atom.sitnam = self.atm_types[ckey].sitnam
                        except AttributeError:
                            pass

        else:
            raise RuntimeError("***Warning No keyword chosen!")

    def fetch_molecules_by_bonds(self, debug=False):
        """
        Form new molecules by checking the bonds (self.bonds) of every atom.
        This function is intended to be a helper function when a topology-file
        is read.
        Atoms are read with their respective ids, which will be converted to
        idx afterwards.
        """
        self.molecules = []  # reset molecules
        tmp_atoms = {}

        if debug is True:
            print("***Info: Combining atoms to molecules by given bond-affiliation.")

        # first we need a dictionary where we can later look up the molecule
        # membership of each atom; initially every atom belongs to None
        for i in self.atoms:
            tmp_atoms[i.atm_id] = None

        # check bonds and assign molecules to connected atoms
        mol_id = 0
        for cbnd in self.bonds:
            atm_1 = cbnd.atm_id1
            atm_2 = cbnd.atm_id2

            # differentiate four cases:
            # 1.) new molecule, maybe part of existing molecule
            if tmp_atoms[atm_1] is None and tmp_atoms[atm_2] is None:
                tmp_atoms[atm_1] = mol_id
                tmp_atoms[atm_2] = mol_id
                self.molecules.append([atm_1, atm_2])
                mol_id += 1
            # 2.) same molecule
            elif tmp_atoms[atm_1] == tmp_atoms[atm_2]:
                if atm_2 not in self.molecules[tmp_atoms[atm_1]]:
                    self.molecules[tmp_atoms[atm_1]].append(atm_2)
            # 3.) same molecule, new atom 1
            elif (tmp_atoms[atm_1] != tmp_atoms[atm_2]) and tmp_atoms[atm_2] is None:
                tmp_atoms[atm_2] = tmp_atoms[atm_1]
                self.molecules[tmp_atoms[atm_1]].append(atm_2)
            # 3.) same molecule, new atom 2
            elif (tmp_atoms[atm_1] != tmp_atoms[atm_2]) and tmp_atoms[atm_1] is None:
                tmp_atoms[atm_1] = tmp_atoms[atm_2]
                self.molecules[tmp_atoms[atm_2]].append(atm_1)
            # 4.) same molecule, currently members of different molecules
            elif tmp_atoms[atm_1] != tmp_atoms[atm_2]:
                # temporary index to overwrite current entry (later)
                tmp = tmp_atoms[atm_2]

                for cur_atm_id in self.molecules[tmp_atoms[atm_2]]:

                    if cur_atm_id not in self.molecules[tmp_atoms[atm_1]]:
                        # only append if necessary
                        self.molecules[tmp_atoms[atm_1]].append(cur_atm_id)

                    # change molecule-ids to be the same as atm_1
                    tmp_atoms[cur_atm_id] = tmp_atoms[atm_1]
                # overwrite entry since no atoms are left; yet the indices must
                # remain intact!
                self.molecules[tmp] = None
            else:
                pass

        # check for atoms that do not form any bonds, e.g. single atom ions
        for catm in tmp_atoms:
            if tmp_atoms[catm] is None:
                self.molecules.append([catm])

        # delete all remaining Nones from self.molecules
        self.molecules = [set(i) for i in self.molecules if i is not None]

    def mols_to_grps(self, debug=False):
        """
        Convert found molecules from fetch_molecules_by_bonds-method to atoms.grp-attribute from the Atoms-Class which
        can be especially helpful for lammps or faster atom-by-molecule search.
        """
        if debug is True:
            print("***Info: Assigning group-attribute of atoms to corresponding molecule indices.")

        for cidx, cmol in enumerate(self.molecules):
            for catm_id in cmol:
                # assign current group-id to current molecule-index
                self.atoms[catm_id].grp_id = cidx

    def same_molecule_as(self, subarray=False, *atm_idxs):
        """
        Get the whole molecule from given atom-indices.
        Input:
            > atm_ids   list or ints; atom-idx (list indices, not ids!) of desired atoms
            > subarray  boolean; output is list of list instead of list of ints;
                        all atoms of each molecule are in one sub-array (flatten array)

        Returns:
            > fetched_molecules list; all atom-idx of found molecules which belong
                                to the same molecule
        """
        # new
        atm_idxs = set(atm_idxs)

        fetched_idxs = []
        for mol_idx, mol in enumerate(self.molecules):
            for atm in mol:
                if atm in atm_idxs:
                    if subarray is False:
                        fetched_idxs.extend(self.molecules[mol_idx])
                    else:
                        fetched_idxs.append(self.molecules[mol_idx])
                    break

        return fetched_idxs

    def add_topology_replicate(self, n, refresh_bonds=False):
        """
        Replicate the given topology n-times (useful when building super cells.)
        Replicates the atoms, bonds, etc. entries.
        """
        print("***Info: Replicating topology")
        atms_old_length = len(self.atoms)
        bnds_old_length = len(self.bonds)
        angs_old_length = len(self.angles)
        dihs_old_length = len(self.dihedrals)
        imps_old_length = len(self.impropers)

        # n-1 since we are adding these to an already existing topology
        for i in xrange(0, n):

            last_atm_id = len(self.atoms)
            # always alter the same first atoms
            for catm in self.atoms[:atms_old_length]:
                # copy the instance, since this is what we want
                catm = copy.copy(catm)
                catm.atm_id += last_atm_id
                self.atoms.append(catm)

            last_bnd_id = len(self.bonds)
            for cbnd in self.bonds[:bnds_old_length]:
                cbnd = copy.copy(cbnd)
                cbnd.bnd_id  += last_bnd_id
                cbnd.atm_id1 += last_atm_id
                cbnd.atm_id2 += last_atm_id
                self.bonds.append(cbnd)

            last_ang_id = len(self.angles)
            for cang in self.angles[:angs_old_length]:
                cang = copy.copy(cang)
                cang.ang_id  += last_ang_id
                cang.atm_id1 += last_atm_id
                cang.atm_id2 += last_atm_id
                cang.atm_id3 += last_atm_id
                self.angles.append(cang)

            last_dih_id = len(self.dihedrals)
            for cdih in self.dihedrals[:dihs_old_length]:
                cdih = copy.copy(cdih)
                cdih.dih_id  += last_dih_id
                cdih.atm_id1 += last_atm_id
                cdih.atm_id2 += last_atm_id
                cdih.atm_id3 += last_atm_id
                cdih.atm_id4 += last_atm_id
                self.dihedrals.append(cdih)

            last_imp_id = len(self.impropers)
            for cimp in self.impropers[:imps_old_length]:
                cimp = copy.copy(cimp)
                cimp.imp_id  += last_imp_id
                cimp.atm_id1 += last_atm_id
                cimp.atm_id2 += last_atm_id
                cimp.atm_id3 += last_atm_id
                cimp.atm_id4 += last_atm_id
                self.impropers.append(cimp)

        if refresh_bonds is True:
            self.fetch_molecules_by_bonds()

    def change_indices(self, incr=1, mode="increase",
                      entries="atm_id, atm_grp_id, atm_key, "):
        """
        Programs, such as lammps, need (why so ever) to have the integers in
        data-files to run from 1-N. Since other programs (like VMD) have starting
        indices with 0, it is easier to convert the internal structure (starting
        with 0) to the one desired.

        Input:
            > mode      str; increase|decrease
        """
        if mode == "increase":
            mod = incr
        elif mode == "decrease":
            mod = -1*incr
        else:
            raise Warning("***Error: Unknown mode '{}'. 'increase' or 'decrease only!'".format(mode))

        for tmp_entry in ["atm_types", "bnd_types", "ang_types", "dih_types", "imp_types"]:
            tmp_dict = {}

            # simplify (saves us repeating the same stuff over and over)
            if tmp_entry == "atm_types":
                cet = self.atm_types
            elif tmp_entry == "bnd_types":
                cet = self.bnd_types
            elif tmp_entry == "ang_types":
                cet = self.ang_types
            elif tmp_entry == "dih_types":
                cet = self.dih_types
            elif tmp_entry == "imp_types":
                cet = self.imp_types
            else:
                pass

            # convert keys
            for cid in cet:
                tmp_dict[cid+mod] = cet[cid]

            if tmp_entry == "atm_types":
                self.atm_types = tmp_dict
            elif tmp_entry == "bnd_types":
                self.bnd_types = tmp_dict
            elif tmp_entry == "ang_types":
                self.ang_types = tmp_dict
            elif tmp_entry == "dih_types":
                self.dih_types = tmp_dict
            else:
                self.imp_types = tmp_dict

        # pair types
        for idx in xrange(len(self.pair_types)):
            self.pair_types[idx].atm_key_i += mod
            #if hasattr(self.pair_types, "atm_key_j"):
            #    self.pair_types[idx].atm_key_j += 1
            self.pair_types[idx].atm_key_j += mod

        # atoms
        for idx in xrange(len(self.atoms)):
            if hasattr(self.atoms[idx], "atm_id") and "atm_id" in entries:
                self.atoms[idx].atm_id  += mod

            if hasattr(self.atoms[idx], "grp_id"):
                self.atoms[idx].grp_id  += mod

            if hasattr(self.atoms[idx], "atm_key"):
                self.atoms[idx].atm_key += mod

        # bonds
        for idx in xrange(len(self.bonds)):
            if hasattr(self.bonds[idx], "bnd_id"):
                self.bonds[idx].bnd_id  += mod
            if hasattr(self.bonds[idx], "bnd_key"):
                self.bonds[idx].bnd_key += mod
            if hasattr(self.bonds[idx], "atm_id1"):
                self.bonds[idx].atm_id1 += mod
            if hasattr(self.bonds[idx], "atm_id2"):
                self.bonds[idx].atm_id2 += mod

        # angles
        for idx in xrange(len(self.angles)):
            if hasattr(self.angles[idx], "ang_id"):
                self.angles[idx].ang_id   += mod
            if hasattr(self.angles[idx], "ang_key"):
                self.angles[idx].ang_key  += mod
            if hasattr(self.angles[idx], "atm_id1"):
                self.angles[idx].atm_id1  += mod
            if hasattr(self.angles[idx], "atm_id2"):
                self.angles[idx].atm_id2  += mod
            if hasattr(self.angles[idx], "atm_id3"):
                self.angles[idx].atm_id3  += mod

        # dihedrals
        for idx in xrange(len(self.dihedrals)):
            if hasattr(self.dihedrals[idx], "dih_id"):
                self.dihedrals[idx].dih_id  += mod
            if hasattr(self.dihedrals[idx], "dih_key"):
                self.dihedrals[idx].dih_key += mod
            if hasattr(self.dihedrals[idx], "atm_id1"):
                self.dihedrals[idx].atm_id1 += mod
            if hasattr(self.dihedrals[idx], "atm_id2"):
                self.dihedrals[idx].atm_id2 += mod
            if hasattr(self.dihedrals[idx], "atm_id3"):
                self.dihedrals[idx].atm_id3 += mod
            if hasattr(self.dihedrals[idx], "atm_id4"):
                self.dihedrals[idx].atm_id4 += mod

        # impropers
        for idx in xrange(len(self.impropers)):
            if hasattr(self.impropers[idx], "imp_id"):
                self.impropers[idx].imp_id  += mod
            if hasattr(self.impropers[idx], "imp_key"):
                self.impropers[idx].imp_key += mod
            if hasattr(self.impropers[idx], "atm_id1"):
                self.impropers[idx].atm_id1 += mod
            if hasattr(self.impropers[idx], "atm_id2"):
                self.impropers[idx].atm_id2 += mod
            if hasattr(self.impropers[idx], "atm_id3"):
                self.impropers[idx].atm_id3 += mod
            if hasattr(self.impropers[idx], "atm_id4"):
                self.impropers[idx].atm_id4 += mod

        # molecules
        for idx in xrange(len(self.molecules)):
            self.molecules[idx] = [i+mod for i in self.molecules[idx]]

        # linked cells
        for frame in xrange(len(self.ts_lnk_cls)):
            for idx_i, i in enumerate(self.ts_lnk_cls[frame].linked_cells):
                for idx_j, j in enumerate(i):
                    for idx_k, k in enumerate(j):
                        #print(idx_i, idx_j, idx_k)
                        self.ts_lnk_cls[frame].linked_cells[idx_i][idx_j][idx_k] = [
                            atm_idx+mod for atm_idx in self.ts_lnk_cls[frame].linked_cells[idx_i][idx_j][idx_k]
                        ]

    def transpose_by_cog(self, frame_id, destination, copy=True):
        """
        Move the current center of geometry so it matches the origin
        at (0/0/0) (all atoms of the same molecule are moved as well).
        Input:
            frame_id        int; id of the frame whose coordinates will be moved
            destination     numpy-array; (1,3)-array with coordinates to move
                            the origin to (inclusive all atoms)
        """
        natoms = len(self.atoms)
        atm_idxs = range(natoms)
        cur_cog = self.get_cog(frame_id, *atm_idxs)
        # 'origin - sys_main_cog' since we want to move the system back
        M_trans_main = cgt.translation_matrix(destination - cur_cog)
        self.mm_atm_coords(0, M_trans_main, copy, *atm_idxs)

    def get_rmsd(self, frame_id_1, frame_id_2, *atm_idxs):
        """
        Calculate the RMSD between two frames.
        Sources: https://pypi.python.org/pypi/rmsd
        """
        P = np.array([self.ts_coords[frame_id_1][i] for i in atm_idxs])
        Q = np.array([self.ts_coords[frame_id_2][i] for i in atm_idxs])
        P -= rmsd.centroid(P)
        Q -= rmsd.centroid(Q)
        pq_rmsd = rmsd.kabsch_rmsd(P, Q)
        return pq_rmsd

    def cut_shape(self, frame_id=-1, inverse_cut=True, *shape_surfaces):
        """
        Cut a certain shape by given shape surfaces. Direction of the normal is important
        since it 'decides' on which side atoms are accepted. Atoms which have
        negative distances (regarding the surface) are outside, atoms with
        positive distances are inside the box. Use the get_plane function from
        the ag_vectalg module to generate suitable planes

        Input:
            > inverse_cut   boolean; decide whether atoms inside (True) or out-
                            side (False) the box are accepted
            > surfaces      tuple; contains all plane-variables (a, b, c, d),
                            mind the sign for each plane since it determines
                            which atom is in- or outside!
        """
        # define some containers for atoms inside and outside the box
        inside_atoms  = []
        outside_atoms = []

        # get atoms that are inside our defined box
        for atm_idx, ccoord in enumerate(self.ts_coords[frame_id]):
            for box_face in shape_surfaces:
                # calculate a simple distance (just interested in the sign)
                distance = agv.get_point_plane_dist(ccoord, *box_face,
                                                    distance_sign_only=True)
                # atom is outside the box -> stop checking the other surfaces
                if distance < 0.0:
                    break

            # atom is inside the box
            if distance >= 0:
                inside_atoms.append(atm_idx)
            else:
                outside_atoms.append(atm_idx)

        if inverse_cut is True:
            # delete atoms which lie outside the box
            self.delete_atoms(*outside_atoms)
        else:
            # delete atoms which are inside the box
            self.delete_atoms(*inside_atoms)

    def find_h_bonds(self, frame_id, distance):
        """
        Find H-Bonds.

        Measure the number of H-Bonds and return a list with atom-indices of all
        atoms involved in the H-Bond, i.e. Acceptor, H-Atom, Donor.
        """

        # create linked cells of non already exist
        if self.ts_lnk_cls != []:
            print("***Info: Creating linked cells.")
            self.create_linked_cells(frame_id)
