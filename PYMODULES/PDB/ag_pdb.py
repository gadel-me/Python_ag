import numpy as np

# import copy
import md_box as mdb
import md_stars as mds
import md_universe as mdu

# import struct
# import ag_lmpdcd_helpers as agldh
# from natsort import natsorted

__version__ = "2019-03-07"


class PdbStuff(mdu.Universe):
    """
    Read file in pdb format. This is very early alpha stage and lots of stuff
    from the pdb format cannot (yet) be read.

    Sources: http://www.wwpdb.org/documentation/file-format-content/format33/sect8.html#CRYST1
    """

    def __init__(self):
        mdu.Universe.__init__(self)

    def read_pdb(self, pdb, overwrite_data=False, debug=False):
        """
        Read a pdb file. Work in Progress!
        """
        with open(pdb, "r") as pdb_file:
            line = pdb_file.readline()
            atm_idx = 0
            bnd_idx = 0
            atm_id_old_new = {}
            all_pdb_coords = []
            tmp_bnds = []

            while line != "":
                if line.startswith("HETATM"):
                    atom_line = line.split()
                    atm_id = int(atom_line[1])
                    sitnam = atom_line[2]
                    coords = np.array([float(i) for i in atom_line[3:6]])
                    atm_id_old_new[atm_id] = atm_idx
                    cur_atm = mds.Atom(atm_id=atm_idx, sitnam=sitnam)
                    self.atoms.append(cur_atm)
                    all_pdb_coords.append(coords)
                    atm_idx += 1
                elif line.startswith("CONECT"):
                    bond_line = line.split()
                    # convert old indices to new ones
                    new_indices = [atm_id_old_new[int(i)] for i in bond_line[1:]]
                    for i, atom_index in enumerate(new_indices):

                        # only covalent bonds are of interest (columns 1-4)
                        if i != 0 and i < 5:
                            id_1 = new_indices[0]
                            id_2 = atom_index

                            # sort by id
                            if id_1 > id_2:
                                id_1, id_2 = atom_index, new_indices[0]

                            if [id_1, id_2] not in tmp_bnds:
                                tmp_bnds.append([id_1, id_2])
                                cbnd = mds.Bond(
                                    bnd_id=bnd_idx, atm_id1=id_1, atm_id2=id_2
                                )
                                bnd_idx += 1
                                self.bonds.append(cbnd)
                elif line.startswith("CRYST1"):
                    a = float(line[7:16])
                    b = float(line[16:25])
                    c = float(line[25:34])
                    alpha = float(line[34:41])
                    beta = float(line[40:47])
                    gamma = float(line[48:55])
                    cbox = mdb.Box(
                        boxtype="lattice",
                        ltc_a=a,
                        ltc_b=b,
                        ltc_c=c,
                        ltc_alpha=alpha,
                        ltc_beta=beta,
                        ltc_gamma=gamma,
                    )
                    self.ts_boxes.append(cbox)
                else:
                    pass

                line = pdb_file.readline()

            # append coordinates
            self.ts_coords.append(all_pdb_coords)

        self.fetch_molecules_by_bonds()
        # self.mols_to_grps()
