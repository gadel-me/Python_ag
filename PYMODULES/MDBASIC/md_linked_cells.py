from __future__ import print_function, division
import time
import math
import numpy as np
import ag_cryst as agc

__version__ = "2019-04-30"


class LinkedCells(object):
    """
    Create linked cells for the atoms of a simulation box.
    """
    def __init__(self,
                 atm_coords,
                 ltc_a, ltc_b, ltc_c,
                 ltc_alpha, ltc_beta, ltc_gamma,
                 coords_type="cartesian"):
        """
        In order to create linked cells, the coordinates of the atoms
        and the box are needed.

        Input:
        """
        self.ltc_a = ltc_a
        self.ltc_b = ltc_b
        self.ltc_c = ltc_c
        self.ltc_alpha = ltc_alpha
        self.ltc_beta  = ltc_beta
        self.ltc_gamma = ltc_gamma
        self.ra = None
        self.rb = None
        self.rc = None

        # convert cartesian-coordinates to fractional
        if coords_type != "fractional":
            atm_coords = self._convert_coords(atm_coords, to_fractional=True)

        self.atm_coords = atm_coords

    def _convert_coords(self, atm_coords, to_fractional=False, to_cartesian=False):
        """
        Convert coordinates to cartesian or fractional coordinates.
        """
        if to_fractional is True:
            # matrix for conversion: cartesian -> fractional
            M_cf = agc.M_cart2fract(self.ltc_a, self.ltc_b, self.ltc_c,
                                    self.ltc_alpha, self.ltc_beta,
                                    self.ltc_gamma)
        elif to_cartesian is True:
            # matrix for conversion: fractional -> cartesian
            M_cf = agc.M_fract2cart(self.ltc_a, self.ltc_b, self.ltc_c,
                                    self.ltc_alpha, self.ltc_beta,
                                    self.ltc_gamma)
        else:
            raise ValueError("fractional or cartesian must be True!")

        # transpose coords for 3x3-matrix multiplication
        ts_coords_T = np.transpose(atm_coords)
        matmul_coords_T = np.matmul(M_cf, ts_coords_T)
        # transpose back
        converted_coords = np.transpose(matmul_coords_T)
        # change coordinates to fractional coordinates
        return converted_coords

    def create_lnk_cells(self, rcut_a=2, rcut_b=2, rcut_c=2, debug=False):
        """
        Divide the cell into ra*rb*rc sub-cells with side-lengths rcut_a|_b|_c

        Data structure:
        e.g. linked_cells[0][0][0]:
            linked_cells[0]        -> sub-cells for a
            linked_cells[0][0]     -> sub-cells for b
            linked_cells[0][0][0]  -> sub-cell for c

        Input:
            > ra, rb, rc    float or int; factors that divide the cell vectors
                            a, b, c into ra, rb, rc sub-cells
                            e.g. a=10, ra=2 -> 5 sub-cells along a
        Returns:
            > linked_cells  list; 3-dimensional array with atom-idx of each sub
                            cell
        """
        if debug is True:
            print("***Linked-Cells Info: Building linked cells.")
            start = time.time()

        # cell division-factors
        self.ra = int(math.ceil(self.ltc_a/rcut_a))
        self.rb = int(math.ceil(self.ltc_b/rcut_b))
        self.rc = int(math.ceil(self.ltc_c/rcut_c))

        if debug is True:
            print("***Linked-Cells Info: Box side lengths had to be adjusted to fit box vectors:\n" +
                  "{:<22s}Side a: {:.3f}, Side b: {:.3f} Side_c: {:.3f}".format(" ",
                                                                                self.ltc_a/self.ra,
                                                                                self.ltc_b/self.rb,
                                                                                self.ltc_c/self.rc))

        avail_max_dist = math.sqrt(rcut_a**2+rcut_b**2+rcut_c**2)

        if debug is True:
            print("***Linked-Cells Info: Max distance between two atoms: {}".format(avail_max_dist))

        linked_cells = []

        # create containers for linked cells
        for sub_a in xrange(self.ra):
            # create first sub row
            clst1 = []

            for sub_b in xrange(self.rb):
                # create second sub row
                clst2 = []

                for sub_c in xrange(self.rc):
                    # create third sub row
                    clst3 = []
                    clst2.append(clst3)

                clst1.append(clst2)

            linked_cells.append(clst1)

        # sub cell fractions
        rca = 1/self.ra
        rcb = 1/self.rb
        rcc = 1/self.rc

        # get all sub lengths (fractional coordinates) in which the cell is divided
        # (for each box vector)
        subcell_a = [i*rca for i in xrange(1, self.ra+1)]  # e.g. 1/3, 2/3, 3/3
        subcell_b = [i*rcb for i in xrange(1, self.rb+1)]
        subcell_c = [i*rcc for i in xrange(1, self.rc+1)]

        # sub cell-index dictionary
        atm_idx_sub_cell = {}

        for cidx in xrange(len(self.atm_coords)):
            # create container for sub cell indices
            atm_idx_sub_cell[cidx] = []

            # wrap coordinates (a(0), b(1) or c(2) > 1 or < 0) back into the cell
            for i in (0, 1, 2):
                # a-coordinate
                if self.atm_coords[cidx][i] > 1:
                    self.atm_coords[cidx][i] -= 1
                elif self.atm_coords[cidx][i] < 0:
                    self.atm_coords[cidx][i] += 1
                else:
                    pass

            # assign atoms to sub cells; csca = current sub cell a
            for csca_idx, csca in enumerate(subcell_a):
                # check if difference between coordinate and current sub cell-
                # vector a is smaller than one sub cell side length rca
                sub_a_distance = csca - self.atm_coords[cidx][0]
                if sub_a_distance >= 0 and sub_a_distance <= rca:
                    break  # stop searching if current atom is inside this cell

            for cscb_idx, cscb in enumerate(subcell_b):
                sub_b_distance = cscb - self.atm_coords[cidx][1]
                if sub_b_distance > 0 and sub_b_distance < rcb:
                    break

            for cscc_idx, cscc in enumerate(subcell_c):
                sub_c_distance = cscc - self.atm_coords[cidx][2]
                if sub_c_distance > 0 and sub_c_distance < rcc:
                    break

            # append atom-index to current box
            linked_cells[csca_idx][cscb_idx][cscc_idx].append(cidx)

        self.linked_cells = linked_cells

        # some verbose stuff
        end = time.time()
        total_cells = self.ra*self.rb*self.rc

        if debug is True:
            print("***Linked-Cells Info: Took {} seconds to build {} linked cells.".format((end - start),
                  total_cells))
