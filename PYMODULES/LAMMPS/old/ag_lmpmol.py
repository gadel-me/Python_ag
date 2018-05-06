"""
Read and write files that follow the lammps molecule standard.
CAVEAT: Coords entry must be the first one to come in order for this module to work.
Sources : http://lammps.sandia.gov/doc/molecule.html
"""
from __future__ import print_function
import numpy as np
import copy
import md_box as mdb
import md_stars as mds
import md_universe as mdu


class LmpMol(mdu.Universe):
    """
    Lammps molecule.
    """
    def __init__(self):
        mdu.Universe.__init__(self)

    def read_lmpmol(self, lmpmol):
        idx_ext_int = {}  # translates between external and internal indices

        with open("lmpmol", "r") as lmpmol_in:
            title = lmpmol_in.next().strip()
            for line in lmpmol:

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
                elif "mass" in line:
                    total_mass = float(line.split()[0])
                elif "com" in line:
                    line = line.split()[0:-1]
                    com = np.array([[float(i) for i in line]])
                elif "inertia" in line:
                    line = line.split()[0:-1]
                    Ixx, Iyy, Izz, Ixy, Ixz, Iyz = [float(i) for i in line]
                elif "Coords" in line:
                    lmpmol_in.next()

                    for atmcnt in xrange(total_atms):
                        line = line.next().split()
                        cid  = int(line[0])
                        coords = np.array([float(i) for i in line[1:4]])

                        self.ts_coords[0].append(coords)
                        self.atoms.append(mds.Atom(atm_id=cid))
                        idx_ext_int[cid] = atmcnt

                elif "Types" in line:
                    lmpmol.next()  # skip empty line

                    for atmcnt in xrange(total_atms):
                        line = lmpmol_in.next().split()
                        cid  = int(line[0])
                        cidx = idx_ext_int[cid]
                        self.atoms[cidx].atm_key = int(line[1])

                elif "Charges" in line:
                    lmpmol.next()  # skip empty line

                    for atmcnt in xrange(total_atms):
                        line = lmpmol_in.next().split()
                        cid  = int(line[0])
                        ccharge = float(line[1])
                        cidx = idx_ext_int[cid]
                        self.atoms[cidx].chge = ccharge

                elif "Diameters" in line:
                    print("Diameters are not implemented (yet).")
                    lmpmol.next()  # skip empty line

                elif "Masses" in line:
                    lmpmol.next()  # skip empty line

                    for atmcnt in xrange(total_atms):
                        line = lmpmol_in.next().split()
                        cid  = int(line[0])
                        cidx = idx_ext_int[cid]
                        cweigh = float(line[1])
                        self.atoms[cidx].weigh = cweigh

                elif "Bonds" in line:
                    lmpmol.next()  # skip empty line

                    for bndcnt in xrange(total_bnds):
                        line = lmpmol_in.next().split()
                        cid, ckey, atm_id1, atm_id2 = [int(i) for i in line[0:4]]
                        atm_idx1 = idx_ext_int[atm_id1]
                        atm_idx2 = idx_ext_int[atm_id2]

                        cbnd = mds.Bond(bnd_id=bndcnt,
                                        bnd_key=ckey,
                                        atm_id1=atm_idx1,
                                        atm_id2=atm_idx2)
                        self.bonds.append(cbnd)

                elif "Angles" in line:
                    lmpmol.next()  # skip empty line

                    for angcnt in xrange(total_angs):
                        line = lmpmol_in.next().split()
                        cid, ckey, atm_id1, atm_id2 = [int(i) for i in line[0:4]]
                        atm_idx1 = idx_ext_int[atm_id1]
                        atm_idx2 = idx_ext_int[atm_id2]

                        cbnd = mds.Bond(bnd_id=bndcnt,
                                        bnd_key=ckey,
                                        atm_id1=atm_idx1,
                                        atm_id2=atm_idx2)
                        self.bonds.append(cbnd)

                elif "Dihedrals" in line:
                    lmpmol.next()  # skip empty line

                elif "Impropers" in line:
                    lmpmol.next()  # skip empty line

                elif "Special Bond Counts" in line:
                    lmpmol.next()  # skip empty line

                elif "Special Bonds" in line:
                    lmpmol.next()  # skip empty line

                elif "Shake Flags" in line:
                    lmpmol.next()  # skip empty line

                elif "Shake Atoms" in line:
                    lmpmol.next()  # skip empty line

                elif "Shake Bond Types" in line:
                    lmpmol.next()  # skip empty line

