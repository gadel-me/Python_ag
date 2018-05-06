from __future__ import print_function
import ag_vectalg as agv
import ag_cryst as agc
#import numpy as np

__version__ = "2017-03-30"


class Box(object):
    """
    Settings for simulation-box.
    """
    def __init__(self, boxtype="lattice",
                 crt_a=None, crt_b=None, crt_c=None,
                 ltc_a=None, ltc_b=None, ltc_c=None,
                 ltc_alpha=None, ltc_beta=None, ltc_gamma=None,
                 lmp_xlo=None, lmp_xhi=None,
                 lmp_ylo=None, lmp_yhi=None,
                 lmp_zlo=None, lmp_zhi=None,
                 lmp_xy=None, lmp_xz=None, lmp_yz=None):
        """
        Get box attributes.
        boxtype: 'cartesian', 'lattice', 'lammps'
        crt_a, crt_b, crt_c: (1,3)-lists of floats/ints; cartesian box vectors
        ltc_a, ltc_b, ltc_c: floats/ints; lattice box vectors
        ltc_alpha, ltc_beta, ltc_gamma: floats/ints; lattice angles
        lmp_xlo, lmp_xhi, lmp_ylo, lmp_yhi, lmp_zlo, lmp_zhi: floats; lammps box size parameters
        lmp_xy, lmp_xz, lmp_yz: floats; lammps tilt factors
        Caveat: ltc_a and crt_a may coexist when a pw input file is read!
                If this is the case, ltc_a is the same as alat or celldm(1)

        Sources:    http://lammps.sandia.gov/doc/Section_howto.html#howto-12
                    https://www.quantum-espresso.org/Doc/INPUT_PW.html#celldm
        """
        self.boxtype = boxtype

        if boxtype == "cartesian":
            self.crt_a = crt_a
            self.crt_b = crt_b
            self.crt_c = crt_c
        elif boxtype == "lattice":
            self.ltc_a = ltc_a
            self.ltc_b = ltc_b
            self.ltc_c = ltc_c
            self.ltc_alpha = ltc_alpha
            self.ltc_beta = ltc_beta
            self.ltc_gamma = ltc_gamma
        elif boxtype == "lammps":
            self.lmp_xlo = lmp_xlo
            self.lmp_xhi = lmp_xhi
            self.lmp_ylo = lmp_ylo
            self.lmp_yhi = lmp_yhi
            self.lmp_zlo = lmp_zlo
            self.lmp_zhi = lmp_zhi
            self.lmp_xy = lmp_xy
            self.lmp_xz = lmp_xz
            self.lmp_yz = lmp_yz
        else:
            raise AttributeError("'boxtype' has to be 'cartesian', " +
                                 "'lattice' or 'lammps'!")

    def box_cart2lat(self):
        """
        Cartesian -> Lattice
        Input: a, b, c = (1,3)-lists
        """
        (self.ltc_a, self.ltc_b, self.ltc_c, self.ltc_alpha, self.ltc_beta,
         self.ltc_gamma) = agc.box_cart2lat(self.crt_a,
                                            self.crt_b,
                                            self.crt_c)
        # change box
        self.boxtype = "lattice"
        del (self.crt_a, self.crt_b, self.crt_c)

    def box_cart2lmp(self):
        """
        Cartesian -> Lammps
        a, b, c = (1,3)-lists
        """
        lx, ly, lz, self.lmp_xy, self.lmp_xz, self.lmp_yz = agc.box_cart2lmp(self.crt_a,
                                                                             self.crt_b,
                                                                             self.crt_c)

        self.lmp_xhi = lx/2
        self.lmp_yhi = ly/2
        self.lmp_zhi = lz/2
        self.lmp_xlo = -self.lmp_xhi
        self.lmp_ylo = -self.lmp_yhi
        self.lmp_zlo = -self.lmp_zhi
        # change box
        self.boxtype = "lammps"
        del (self.crt_a, self.crt_b, self.crt_c)

    def box_lat2cart(self):
        """
        Lattice -> Cartesian
         alpha, beta, gamma = floats (radians)
        """
        self.crt_a, self.crt_b, self.crt_c = agc.box_lat2cart(self.ltc_a,
                                                              self.ltc_b,
                                                              self.ltc_c,
                                                              self.ltc_alpha,
                                                              self.ltc_beta,
                                                              self.ltc_gamma)
        # change box
        self.boxtype = "cartesian"
        del (self.ltc_a, self.ltc_b, self.ltc_c, self.ltc_alpha, self.ltc_beta,
             self.ltc_gamma)

    def box_lat2lmp(self, triclinic=True):
        """
        Lattice -> Lammps

        Input:
            a, b, c              floats; lattice vectors a, b and c
            alpha, beta, gamma   floats (radians); lattice angles alpha, beta, gamma
            triclinic            boolean; always delete (mark as None)
                                 xy, xz and yz in order to make cell orthogonal
        """
        lx, ly, lz, self.lmp_xy, self.lmp_xz, self.lmp_yz = agc.box_lat2lmp(self.ltc_a,
                                                                            self.ltc_b,
                                                                            self.ltc_c,
                                                                            self.ltc_alpha,
                                                                            self.ltc_beta,
                                                                            self.ltc_gamma)
        self.lmp_xhi = lx/2
        self.lmp_yhi = ly/2
        self.lmp_zhi = lz/2
        self.lmp_xlo = -self.lmp_xhi
        self.lmp_ylo = -self.lmp_yhi
        self.lmp_zlo = -self.lmp_zhi
        # change box
        self.boxtype = "lammps"

        # define xy, xz, yz as None (needed for writing lammps-data routine)
        if triclinic is False:
            self.lmp_xy = self.lmp_xz = self.lmp_yz = None

        del (self.ltc_a, self.ltc_b, self.ltc_c, self.ltc_alpha, self.ltc_beta,
             self.ltc_gamma)

    def box_lmp2cart(self):
        """
        Lammps -> Cartesian
        Returns: (1,3)-list with cartesian coordinates
        """
        lx = self.lmp_xhi - self.lmp_xlo
        ly = self.lmp_yhi - self.lmp_ylo
        lz = self.lmp_zhi - self.lmp_zlo

        if self.lmp_xy is None:
            self.lmp_xy = 0
        if self.lmp_xz is None:
            self.lmp_xz = 0
        if self.lmp_yz is None:
            self.lmp_yz = 0

        self.crt_a, self.crt_b, self.crt_c = agc.box_lmp2cart(lx, ly, lz,
                                                              self.lmp_xy,
                                                              self.lmp_xz,
                                                              self.lmp_yz)
        # change box
        self.boxtype = "cartesian"
        del (self.lmp_xhi, self.lmp_yhi, self.lmp_zhi, self.lmp_xlo, self.lmp_ylo,
             self.lmp_zlo, self.lmp_xy, self.lmp_xz, self.lmp_yz)

    def box_lmp2lat(self):
        """
        Lammps -> Lattice
        """
        lx = self.lmp_xhi - self.lmp_xlo
        ly = self.lmp_yhi - self.lmp_ylo
        lz = self.lmp_zhi - self.lmp_zlo

        if self.lmp_xy is None:
            self.lmp_xy = 0
        if self.lmp_xz is None:
            self.lmp_xz = 0
        if self.lmp_yz is None:
            self.lmp_yz = 0

        (self.ltc_a, self.ltc_b, self.ltc_c, self.ltc_alpha, self.ltc_beta,
         self.ltc_gamma) = agc.box_lmp2lat(lx, ly, lz,
                                           self.lmp_xy, self.lmp_xz, self.lmp_yz)
        # change box
        self.boxtype = "lattice"
        del (self.lmp_xhi, self.lmp_yhi, self.lmp_zhi, self.lmp_xlo, self.lmp_ylo,
             self.lmp_zlo, self.lmp_xy, self.lmp_xz, self.lmp_yz)

    def get_center(self):
        """
        Works for parallelepiped shaped cubes.
        """
        if self.boxtype == "cartesian":
            self.center = [
                agv.vs_mult(self.crt_a, 0.5),
                agv.vs_mult(self.crt_b, 0.5),
                agv.vs_mult(self.crt_c, 0.5)
            ]
        elif self.boxtype == "lattice":
            self.center = [0.5, 0.5, 0.5]
        elif self.boxtype == "lammps":
            self.center = [
                (abs(self.lmp_xlo)+abs(self.lmp_xhi))*0.5,
                (abs(self.lmp_ylo)+abs(self.lmp_yhi))*0.5,
                (abs(self.lmp_zlo)+abs(self.lmp_zhi))*0.5
            ]
        else:
            pass
