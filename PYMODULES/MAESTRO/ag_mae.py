from __future__ import print_function, division
import numpy as np
import re
import md_stars as mds
import md_universe as mdu

__version__ = "2017-07-26"


class MaestroStuff(mdu.Universe):
    """
    Read Schroedinger Maestro's files.
    """
    def __init__(self):
        """
        Load universe first.
        """
        mdu.Universe.__init__(self)

    def read_mae(self, maefile, overwrite_data=False):
        """
        Read Schroedinger Maestro's mae-files.
        """
        atm_id_old_new = {}

        with open(maefile, "r") as mae_in:
            for line in mae_in:

                # find header for atoms-section
                if "m_atom" in line:
                    cframe = []
                    num_atms = int(re.findall(r'\d+', line)[0])

                    # skip lines until coordinates section is reached
                    while ":::" not in line:
                        line = mae_in.next()

                    for iid in xrange(num_atms):
                        catm = mae_in.next().split()
                        atm_id_old_new[int(catm[0])] = iid
                        csitnam = catm[-1]
                        ccoords = np.array([float(i) for i in catm[2:5]])

                        # check if an instance of Atom with index iid already exists;
                        # overwrite info if it does or create a new one if it does not
                        try:
                            self.atoms[iid]

                            # overwrite data
                            if overwrite_data is True:
                                self.atoms[iid].atm_id = iid
                                self.atoms[iid].sitnam = csitnam
                            # complement data
                            else:

                                if not hasattr(self.atoms[iid], "atm_id"):
                                    self.atoms[iid].atm_id = iid

                                if not hasattr(self.atoms[iid], "sitnam"):
                                    self.atoms[iid].sitnam = csitnam

                        except IndexError:
                                catm = mds.Atom(atm_id=iid,
                                                sitnam=csitnam)
                                self.atoms.append(catm)

                        cframe.append(ccoords)
                    # append current frame
                    self.ts_coords.append(cframe)

                elif "m_bond" in line:
                    num_bnds = int(re.findall(r'\d+', line)[0])

                    # skip lines until coordinates section is reached
                    while ":::" not in line:
                        line = mae_in.next()

                    for iid in xrange(num_bnds):
                        cbnd = mae_in.next().split()

                        # store section
                        atm_1     = atm_id_old_new[int(cbnd[1])]
                        atm_2     = atm_id_old_new[int(cbnd[2])]
                        bnd_order = int(cbnd[3])

                        # check if an instance of Bond with index iid already exists;
                        # overwrite info if it does or create a new one if it does not
                        try:
                            self.bonds[iid]

                            # overwrite data
                            if overwrite_data is True:
                                self.bonds[iid].bnd_id = iid
                                self.bonds[iid].bnd_order = bnd_order
                            # complement data
                            else:

                                if not hasattr(self.bonds[iid], "bnd_id"):
                                    self.bonds[iid].bnd_id = iid

                                if not hasattr(self.bonds[iid], "atm_id1"):
                                    self.bonds[iid].atm_id1 = atm_1

                                if not hasattr(self.bonds[iid], "atm_id2"):
                                    self.bonds[iid].atm_id2 = atm_2

                                if not hasattr(self.bonds[iid], "bnd_order"):
                                    self.bonds[iid].bnd_order = bnd_order

                        except IndexError:
                                cbnd = mds.Bond(bnd_id=iid,
                                                atm_id1=atm_1,
                                                atm_id2=atm_2,
                                                bnd_order=bnd_order)
                                self.bonds.append(cbnd)
                else:
                    pass
