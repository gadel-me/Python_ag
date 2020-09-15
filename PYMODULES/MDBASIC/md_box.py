# import numpy as np
import scipy.constants as sc
import ag_vectalg as agv
import ag_cryst as agc

__version__ = "2017-03-30"

BOHR_ANGSTROM = sc.value("Bohr radius") / sc.angstrom
ANGSTROM_BOHR = sc.angstrom / sc.value("Bohr radius")


class Box(object):
    """
    Settings for simulation-box.
    """

    def __init__(
        self,
        boxtype="lattice",
        unit="angstrom",
        crt_a=None,
        crt_b=None,
        crt_c=None,
        ltc_a=None,
        ltc_b=None,
        ltc_c=None,
        ltc_alpha=None,
        ltc_beta=None,
        ltc_gamma=None,
        lmp_xlo=None,
        lmp_xhi=None,
        lmp_ylo=None,
        lmp_yhi=None,
        lmp_zlo=None,
        lmp_zhi=None,
        lmp_xy=None,
        lmp_xz=None,
        lmp_yz=None,
        alat=None,
    ):
        """
        Get box attributes.
        boxtype: 'cartesian', 'lattice', 'lammps'
        crt_a, crt_b, crt_c: (1,3)-lists of floats/ints; cartesian box vectors
        ltc_a, ltc_b, ltc_c: floats/ints; lattice box vectors
        ltc_alpha, ltc_beta, ltc_gamma: floats/ints; lattice angles
        lmp_xlo, lmp_xhi, lmp_ylo, lmp_yhi, lmp_zlo, lmp_zhi: floats; lammps box size parameters
        lmp_xy, lmp_xz, lmp_yz: floats; lammps tilt factors
        unit    str; angstrom | bohr | alat
                    alat: lattice vector a in bohr (only used in pw), see pw.x input description

        #TODO boxtype should always be given or an error for the user has to
        #TODO appear, that states that a boxtype must be defined!

        Sources:    http://lammps.sandia.gov/doc/Section_howto.html#howto-12
                    https://www.quantum-espresso.org/Doc/INPUT_PW.html#celldm
        """
        self.boxtype = boxtype
        self.unit = unit
        self.volume = None

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
        elif boxtype == "alat":
            self.crt_a = crt_a
            self.crt_b = crt_b
            self.crt_c = crt_c
            self.alat = alat
        else:
            raise AttributeError(
                "'boxtype' has to be 'cartesian', 'alat', " + "'lattice' or 'lammps'!"
            )

    def box_cart2lat(self):
        """
        Cartesian -> Lattice
        Input: a, b, c = (1,3)-lists
        """
        (
            self.ltc_a,
            self.ltc_b,
            self.ltc_c,
            self.ltc_alpha,
            self.ltc_beta,
            self.ltc_gamma,
        ) = agc.box_cart2lat(self.crt_a, self.crt_b, self.crt_c)
        # change box
        self.boxtype = "lattice"
        del (self.crt_a, self.crt_b, self.crt_c)

    def box_cart2lmp(self):
        """
        Cartesian -> Lammps
        a, b, c = (1,3)-lists
        """
        lx, ly, lz, self.lmp_xy, self.lmp_xz, self.lmp_yz = agc.box_cart2lmp(
            self.crt_a, self.crt_b, self.crt_c
        )

        self.lmp_xhi = lx / 2
        self.lmp_yhi = ly / 2
        self.lmp_zhi = lz / 2
        self.lmp_xlo = None if self.lmp_xhi is None else -1 * self.lmp_xhi
        self.lmp_ylo = None if self.lmp_yhi is None else -1 * self.lmp_yhi
        self.lmp_zlo = None if self.lmp_zhi is None else -1 * self.lmp_zhi
        # change box1 *
        self.boxtype = "lammps"
        del (self.crt_a, self.crt_b, self.crt_c)

    def box_lat2cart(self):
        """
        Lattice -> Cartesian
         alpha, beta, gamma = floats (radians)
        """
        self.crt_a, self.crt_b, self.crt_c = agc.box_lat2cart(
            self.ltc_a,
            self.ltc_b,
            self.ltc_c,
            self.ltc_alpha,
            self.ltc_beta,
            self.ltc_gamma,
        )
        # change box
        self.boxtype = "cartesian"
        del (
            self.ltc_a,
            self.ltc_b,
            self.ltc_c,
            self.ltc_alpha,
            self.ltc_beta,
            self.ltc_gamma,
        )

    def box_lat2lmp(self, triclinic=True):
        """
        Lattice -> Lammps

        Input:
            a, b, c              floats; lattice vectors a, b and c
            alpha, beta, gamma   floats (radians); lattice angles alpha, beta, gamma
            triclinic            boolean; always delete (mark as None)
                                 xy, xz and yz in order to make cell orthogonal
        """
        lx, ly, lz, self.lmp_xy, self.lmp_xz, self.lmp_yz = agc.box_lat2lmp(
            self.ltc_a,
            self.ltc_b,
            self.ltc_c,
            self.ltc_alpha,
            self.ltc_beta,
            self.ltc_gamma,
        )

        self.lmp_xhi = lx / 2
        self.lmp_yhi = ly / 2
        self.lmp_zhi = lz / 2
        self.lmp_xlo = (
            None if self.lmp_xhi is None else -1 * self.lmp_xhi
        )  # // -self.lmp_xhi
        self.lmp_ylo = (
            None if self.lmp_yhi is None else -1 * self.lmp_yhi
        )  # // -self.lmp_yhi
        self.lmp_zlo = (
            None if self.lmp_zhi is None else -1 * self.lmp_zhi
        )  # // -self.lmp_zhi
        # change box
        self.boxtype = "lammps"

        # define xy, xz, yz as None (needed for writing lammps-data routine)
        if triclinic is False:
            self.lmp_xy = self.lmp_xz = self.lmp_yz = None

        del (
            self.ltc_a,
            self.ltc_b,
            self.ltc_c,
            self.ltc_alpha,
            self.ltc_beta,
            self.ltc_gamma,
        )

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

        self.crt_a, self.crt_b, self.crt_c = agc.box_lmp2cart(
            lx, ly, lz, self.lmp_xy, self.lmp_xz, self.lmp_yz
        )
        # change box
        self.boxtype = "cartesian"
        del (
            self.lmp_xhi,
            self.lmp_yhi,
            self.lmp_zhi,
            self.lmp_xlo,
            self.lmp_ylo,
            self.lmp_zlo,
            self.lmp_xy,
            self.lmp_xz,
            self.lmp_yz,
        )

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

        (
            self.ltc_a,
            self.ltc_b,
            self.ltc_c,
            self.ltc_alpha,
            self.ltc_beta,
            self.ltc_gamma,
        ) = agc.box_lmp2lat(lx, ly, lz, self.lmp_xy, self.lmp_xz, self.lmp_yz)
        # change box
        self.boxtype = "lattice"
        del (
            self.lmp_xhi,
            self.lmp_yhi,
            self.lmp_zhi,
            self.lmp_xlo,
            self.lmp_ylo,
            self.lmp_zlo,
            self.lmp_xy,
            self.lmp_xz,
            self.lmp_yz,
        )

    def box_cart2alat(self):
        assert self.boxtype == "cartesian"
        self.crt_a, self.crt_b, self.crt_c, self.alat = agc.box_cart2alat(
            self.crt_a, self.crt_b, self.crt_c
        )
        self.boxtype = "alat"

    def box_lmp2alat(self):
        """box_lmp2alat [summary]
        Lammps -> alat
        """
        assert self.boxtype == "lammps"

        lx = self.lmp_xhi - self.lmp_xlo
        ly = self.lmp_yhi - self.lmp_ylo
        lz = self.lmp_zhi - self.lmp_zlo

        if self.lmp_xy is None:
            self.lmp_xy = 0
        if self.lmp_xz is None:
            self.lmp_xz = 0
        if self.lmp_yz is None:
            self.lmp_yz = 0

        self.crt_a, self.crt_b, self.crt_c, self.alat = agc.box_lmp2alat(
            lx, ly, lz, self.lmp_xy, self.lmp_xz, self.lmp_yz
        )
        self.boxtype = "alat"

    def get_center(self):
        """
        Works for parallelepiped shaped cubes.
        """
        if self.boxtype == "cartesian":
            self.center = [
                agv.vs_mult(self.crt_a, 0.5),
                agv.vs_mult(self.crt_b, 0.5),
                agv.vs_mult(self.crt_c, 0.5),
            ]
        elif self.boxtype == "lattice":
            self.center = [0.5, 0.5, 0.5]
        elif self.boxtype == "lammps":
            self.center = [
                (abs(self.lmp_xlo) + abs(self.lmp_xhi)) * 0.5,
                (abs(self.lmp_ylo) + abs(self.lmp_yhi)) * 0.5,
                (abs(self.lmp_zlo) + abs(self.lmp_zhi)) * 0.5,
            ]
        else:
            pass

    def angstrom2bohr(self):
        """
        Change the lengths of each box vector from angstrom to bohr.
        """
        # first check that units are not bohr already
        if self.unit != "bohr":
            if self.boxtype == "lattice":
                self.ltc_a *= ANGSTROM_BOHR
                self.ltc_b *= ANGSTROM_BOHR
                self.ltc_c *= ANGSTROM_BOHR
            elif self.boxtype == "cartesian":
                self.crt_a = [i * ANGSTROM_BOHR for i in self.crt_a]
                self.crt_b = [i * ANGSTROM_BOHR for i in self.crt_b]
                self.crt_c = [i * ANGSTROM_BOHR for i in self.crt_c]
            else:
                raise Warning("Conversion of lammps boxes currently not supported!")

        self.unit = "bohr"

    def bohr2angstrom(self):
        """
        Change the lengths of each box vector from bohr to angstrom
        """
        # first check that units are not bohr already
        if self.unit != "bohr":
            if self.boxtype == "lattice":
                self.ltc_a *= BOHR_ANGSTROM
                self.ltc_b *= BOHR_ANGSTROM
                self.ltc_c *= BOHR_ANGSTROM
            elif self.boxtype == "cartesian":
                self.crt_a = [i * BOHR_ANGSTROM for i in self.crt_a]
                self.crt_b = [i * BOHR_ANGSTROM for i in self.crt_b]
                self.crt_c = [i * BOHR_ANGSTROM for i in self.crt_c]
            else:
                raise Warning("Conversion of lammps boxes currently not supported!")

        self.unit = "angstrom"

    def alat2angstrom(self, alat):
        """
        Convert unit 'alat' which is in bohr (see pw.x input description) to unit angstrom.
        """
        # just check that the current unit is alat
        if self.unit == "alat":
            # check if the box type is cartesian or it will not be converted
            if self.boxtype == "cartesian":
                self.crt_a = [i * BOHR_ANGSTROM * alat for i in self.crt_a]
                self.crt_b = [i * BOHR_ANGSTROM * alat for i in self.crt_b]
                self.crt_c = [i * BOHR_ANGSTROM * alat for i in self.crt_c]
            else:
                print(
                    "Wrong unit with wrong boxtype ({}/{})!".format(
                        self.unit, self.boxtype
                    )
                )
        else:
            print("Wrong unit ({})!".format(self.unit))

        self.unit = "angstrom"

    def calc_volume(self):
        """
        Calculate the volume of the cell.
        """
        self.volume = agc.box_lat_volume(
            self.ltc_a,
            self.ltc_b,
            self.ltc_c,
            self.ltc_alpha,
            self.ltc_beta,
            self.ltc_gamma,
        )
