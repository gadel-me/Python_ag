import math
import md_stars_helper as mdsh
import md_elements as mde

__version__ = "2017-04-10"


class IterMixin(object):
    """
    Class to return the attributes with their corresponding values of the
    desired object.
    """

    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value


class Atom(IterMixin):
    """
    Class for Atoms as Balls and for ForceField-Definitions of AtomTypes
    Van der Waals radii definitions.
    Amber-FF: Energy in kcal/mol!
    Sources:    https://www.researchgate.net/post/Where_are_the_charge_sigma_epsilon_parameters_in_GAFF
                https://docs.scipy.org/doc/scipy/reference/constants.html#scipy.constants.physical_constants
    Currently implemented unit-conversion: kcal/mol, kj/mol, eV
    Remember: r0 is not sigma; r0 = 2**(1/6)*sigma (using Lennard-Jones Potential)
    #TODO CHANGE FROZEN TO A LIST WHERE THE FIRST MEMBER IS THE PROGRAM
    #TODO IT HAS THE INFO FROM AND THEN X Y AND Z AS PARAMETERS ON WHAT TO FREEZE
    #TODO FOR GAUSSIAN X Y AND Z SHALL BE THE SAME NUMBER
    #TODO A CHECK SHOULD BE DONE TO ASSERT THAT THE INFO FITS THE PROGRAM IT WAS FROM
    """

    def __init__(
        self,
        atm_id=None,
        atm_key=None,
        sitnam=None,
        weigh=None,
        grp_id=None,
        chge=None,
        coords=None,
        res=None,
        nrept=None,
        ifrz=None,
        ifrz_x=None,
        ifrz_y=None,
        ifrz_z=None,
        igrp=None,
        r0=None,
        sigma=None,
        epsilon=None,
        energy_unit=None,
        comment=None,
        force=None,
        velocity=None,
        pseudopotential=None,
        ucell_info=None
    ):
        if atm_id is not None:
            self.atm_id = atm_id

        if atm_key is not None:
            self.atm_key = atm_key

        if sitnam is not None:
            self.sitnam = sitnam

        if weigh is not None:
            self.weigh = weigh

        if grp_id is not None:
            self.grp_id = grp_id

        if chge is not None:
            self.chge = chge

        if coords is not None:
            self.coords = coords

        if res is not None:
            self.res = res

        if nrept is not None:
            self.nrept = nrept

        if ifrz is not None:
            self.ifrz = ifrz

        # pw distinguishes freezing
        if ifrz_x is not None:
            self.ifrz_x = ifrz_x

        if ifrz_y is not None:
            self.ifrz_y = ifrz_y

        if ifrz_z is not None:
            self.ifrz_z = ifrz_z

        if igrp is not None:
            self.igrp = igrp

        if r0 is not None:
            self.r0 = r0
            self.sigma = r0 / (2 ** (1 / 6))

        if sigma is not None:
            self.sigma = sigma

        if epsilon is not None:
            self.epsilon = epsilon

        if energy_unit is not None:
            self.energy_unit = energy_unit

        if comment is not None:
            self.comment = comment

        if force is not None:
            self.force = force

        if velocity is not None:
            self.velocity = velocity

        if pseudopotential is not None:
            self.pseudopotential = pseudopotential

        # ids of the unit cell the atom is in, starting with [0, 0, 0]
        #  e.g. 1 2 5 means the atom is in unit cell repetition
        # a * 1 (2nd cell), b * 2 (3rd), c * 5 (6th) counting from
        # the first one (0 0 0)
        # only apllies to super cells consisting of unit cells!
        if ucell_info is not None:
            self.ucell_info = ucell_info

    def convert_energy_unit(self, unit_out):
        """
        Convert the energy unit to the one desired.
        """
        self.epsilon = mdsh.convert_energy_unit(
            self.epsilon, self.energy_unit, unit_out
        )
        self.energy_unit = unit_out

    def mix_ij(self, sigma_j, epsilon_j, mix="arithmetic"):
        """
        Calculate sigma ij from the sigmas of the two atoms i and j.
        Sources:    http://lammps.sandia.gov/doc/pair_modify.html
        """
        epsilon_ij = math.sqrt(self.epsilon * epsilon_j)

        if mix == "arithmetic":
            sigma_ij = (self.sigma + sigma_j) / 2
        elif mix == "geometric":
            sigma_ij = math.sqrt(self.sigma * sigma_j)
        else:
            raise RuntimeError("'mix' has to be 'arithmetic' or 'geometric'!")

        return (sigma_ij, epsilon_ij)

    def calc_weigh(self, overwrite=False):
        """
        Calculate the element mass from the elements name.
        """
        if not hasattr(self, "weigh"):
            self.weigh = None

        if self.weigh is None or overwrite is True:
            try:
                self.weigh = mde.element_mass[self.sitnam]
            except KeyError:
                # check if sitnam is the atomic number
                atomic_number = int(self.sitnam)
                self.weigh = mde.atomic_number_mass[atomic_number]


class Bond(IterMixin):
    """
    Sources:    http://cbio.bmt.tue.nl/pumma/index.php/Theory/Potentials

    prm1 = force constant
    prm2 = optimal distance (minimum distance, r0)
    """

    def __init__(
        self,
        bnd_id=None,
        bnd_key=None,
        atm_id1=None,
        atm_id2=None,
        atm_sitnam_1=None,
        atm_sitnam_2=None,
        prm1=None,
        prm2=None,
        prm3=None,
        prm4=None,
        energy_unit=None,
        bnd_order=None,
        comment=None,
    ):
        """
        bnd_order   float; bond-order; single (1.0)-, double(2.0)-, triple(2.0)-bond
        """
        if bnd_id is not None:
            self.bnd_id = bnd_id
        if bnd_key is not None:
            self.bnd_key = bnd_key
        if atm_id1 is not None:
            self.atm_id1 = atm_id1
        if atm_id2 is not None:
            self.atm_id2 = atm_id2
        if atm_sitnam_1 is not None:
            self.atm_sitnam_1 = atm_sitnam_1
        if atm_sitnam_2 is not None:
            self.atm_sitnam_2 = atm_sitnam_2
        if prm1 is not None:
            self.prm1 = prm1
        if prm2 is not None:
            self.prm2 = prm2
        if prm3 is not None:
            self.prm3 = prm3
        if prm4 is not None:
            self.prm4 = prm4
        if energy_unit is not None:
            self.energy_unit = energy_unit
        if bnd_order is not None:
            self.bnd_order = bnd_order

        self.comment = comment

    def convert_energy_unit(self, unit_out):
        self.prm1 = mdsh.convert_energy_unit(self.prm1, self.energy_unit, unit_out)
        self.energy_unit = unit_out

    def check_bnd_type(self):
        """
        Check if something is odd with the force-field parameters (currently
        only eV-units are supported).
        """
        warning_message = (
            "***Warning: Odd bond-type-parameters found. Check if parameter is o.k."
        )

        if self.prm1 > 30 or self.prm2 > 2:
            print(warning_message)


class Angle(IterMixin):
    """
    Help on class Angle in module md_molecule:
    Sources:    http://lammps.sandia.gov/doc/angle_harmonic.html
                http://cbio.bmt.tue.nl/pumma/index.php/Theory/Potentials

    class Angle(__builtin__.object)
        Angle(ang_id=None, ang_key=None, atm_id1=None, atm_id2=None,
              atm_id3=None, prm1=None, prm2=None, prm3=None, prm4=None)

    prm1 = force constant
    prm2 = minimal angle (optimal angle)
    """

    def __init__(
        self,
        ang_id=None,
        ang_key=None,
        atm_id1=None,
        atm_id2=None,
        atm_id3=None,
        prm1=None,
        prm2=None,
        prm3=None,
        prm4=None,
        energy_unit=None,
        angle_unit=None,
        comment=None,
    ):
        if ang_id is not None:
            self.ang_id = ang_id
        if ang_key is not None:
            self.ang_key = ang_key
        if atm_id1 is not None:
            self.atm_id1 = atm_id1
        if atm_id2 is not None:
            self.atm_id2 = atm_id2
        if atm_id3 is not None:
            self.atm_id3 = atm_id3
        if prm1 is not None:
            self.prm1 = prm1
        if prm2 is not None:
            self.prm2 = prm2
        if prm3 is not None:
            self.prm3 = prm3
        if prm4 is not None:
            self.prm4 = prm4
        if energy_unit is not None:
            self.energy_unit = energy_unit
        if angle_unit is not None:
            self.angle_unit = angle_unit

        self.comment = comment

    def convert_energy_unit(self, unit_out):
        self.prm1 = mdsh.convert_energy_unit(self.prm1, self.energy_unit, unit_out)
        self.energy_unit = unit_out

    def convert_angle_unit(self, unit_out):
        """
        unit_out:   str; 'deg'|'rad'
        """
        self.prm2 = mdsh.convert_angle_unit(self.prm2, unit_out)
        self.angle_unit = unit_out

    def check_ang_type(self):
        """
        Check if something is odd with the force-field parameters (currently
        only eV-units are supported).
        """
        warning_message = (
            "***Warning: Odd angle-type-parameter found. Check if parameter is o.k."
        )

        if self.prm1 > 10:
            print(warning_message)


class Dihedral(IterMixin):
    """
    Help on class Dihedral in module md_molecule:
    Sources:    http://lammps.sandia.gov/doc/dihedral_style.html
                http://lammps.sandia.gov/doc/dihedral_charmm.html
                http://cbio.bmt.tue.nl/pumma/index.php/Theory/Potentials

    prm_k:  dihedral force constant
    prm_n:  dihedral periodicity
    prm_d:  dihedral phase
    """

    def __init__(
        self,
        dih_id=None,
        dih_key=None,
        atm_id1=None,
        atm_id2=None,
        atm_id3=None,
        atm_id4=None,
        prm_k=None,
        prm_n=None,
        prm_d=None,
        weigh_factor=None,
        elec_inter_1_4_scale=None,
        vdw_inter_1_4_scale=None,
        energy_unit=None,
        angle_unit=None,
        comment=None,
    ):
        if dih_id is not None:
            self.dih_id = dih_id
        if dih_key is not None:
            self.dih_key = dih_key
        if atm_id1 is not None:
            self.atm_id1 = atm_id1
        if atm_id2 is not None:
            self.atm_id2 = atm_id2
        if atm_id3 is not None:
            self.atm_id3 = atm_id3
        if atm_id4 is not None:
            self.atm_id4 = atm_id4
        if prm_k is not None:
            self.prm_k = prm_k
        if prm_n is not None:
            self.prm_n = prm_n
        if prm_d is not None:
            self.prm_d = prm_d
        if weigh_factor is not None:
            self.weigh_factor = weigh_factor
        if elec_inter_1_4_scale is not None:
            self.elec_inter_1_4_scale = elec_inter_1_4_scale
        if vdw_inter_1_4_scale is not None:
            self.vdw_inter_1_4_scale = vdw_inter_1_4_scale
        if energy_unit is not None:
            self.energy_unit = energy_unit
        if angle_unit is not None:
            self.angle_unit = angle_unit

        self.comment = comment

    def convert_energy_unit(self, unit_out):
        self.prm_k = mdsh.convert_energy_unit(self.prm_k, self.energy_unit, unit_out)
        self.energy_unit = unit_out

    def convert_angle_unit(self, unit_out):
        """
        unit_out:   str; 'deg'|'rad'
        """
        self.prm_d = mdsh.convert_angle_unit(self.prm_d, unit_out)
        self.angle_unit = unit_out

    def check_dih_type(self):
        """
        Check if something is odd with the force-field parameters (currently
        only eV-units are supported).
        """
        warning_message = (
            "***Warning: Odd dihedral-type-parameter found. Check if parameter is o.k."
        )

        if self.prm_k > 1:
            print(warning_message)

        if int(self.prm_d) not in (180, 0):
            print(warning_message)

    def create_lmp_dih_style(self, dih_style="charmm"):
        """
        Create the lammps 'dihedral_coeff' string.

        Sources:    https://lammps.sandia.gov/doc/dihedral_charmm.html
                    https://lammps.sandia.gov/doc/dihedral_coeff.html

        Parameters
        ----------
        dih_style : str
            dihedral style as in the lammps manual (currently only 'charmm' is supported)

        Returns
        -------
        dih_coeff : str
            dihedral-coeff-string as in the lammps manual
        """
        if self.prm_k is None or self.prm_d is None or self.prm_n is None:
            raise Warning("Force constant k, phase shift d or periodicity n missing!")

        return "dihedral_coeff {} "


class Improper(IterMixin):
    """
    Help on class Improper in module md_molecule:

    Sources (LAMMPS):   http://lammps.sandia.gov/doc/improper_style.html
                        http://lammps.sandia.gov/doc/improper_cvff.html
                        http://cbio.bmt.tue.nl/pumma/index.php/Theory/Potentials
    Sources (DLPOLY):   DLPOLY-MANUAL, p.

    Only cvff-style at the moment:
    prm_k:   improper force constant
    prm_n:   periodicity; (0,1,2,3,4,6)
    prm_d:   -1 or +1 (= cos(0), cos(pi))
    angle_unit: 'rad'|deg
    """

    def __init__(
        self,
        imp_id=None,
        imp_key=None,
        atm_id1=None,
        atm_id2=None,
        atm_id3=None,
        atm_id4=None,
        prm_k=None,
        prm_d=None,
        prm_n=None,
        energy_unit=None,
        angle_unit=None,
        comment=None,
    ):
        if imp_id is not None:
            self.imp_id = imp_id
        if imp_key is not None:
            self.imp_key = imp_key
        if atm_id1 is not None:
            self.atm_id1 = atm_id1
        if atm_id2 is not None:
            self.atm_id2 = atm_id2
        if atm_id3 is not None:
            self.atm_id3 = atm_id3
        if atm_id4 is not None:
            self.atm_id4 = atm_id4
        if prm_k is not None:
            self.prm_k = prm_k
        if prm_n is not None:
            self.prm_n = prm_n
        if prm_d is not None:
            self.prm_d = prm_d
        if energy_unit is not None:
            self.energy_unit = energy_unit
        if angle_unit is not None:
            self.angle_unit = angle_unit

        self.comment = comment

    def convert_energy_unit(self, unit_out):
        self.prm_k = mdsh.convert_energy_unit(self.prm_k, self.energy_unit, unit_out)
        self.energy_unit = unit_out

    def cvff_prm_d(self):
        """
        unit_out:   str; 'deg'|'rad'
        """
        if self.angle_unit != "rad":
            self.prm_d = math.radians(self.prm_d)

        # get cosine
        tmp_prm_d = math.cos(self.prm_d)

        # round to next negative value if cos < 0 (e.g. -0.9999)
        if tmp_prm_d < 0:
            # print("Smaller")
            tmp_prm_d = math.floor(tmp_prm_d)
        # round to next positive int if cos > 0 (e.g. 0.99999)
        else:
            # print("bigger")
            tmp_prm_d = math.ceil(tmp_prm_d)

        self.prm_d = int(tmp_prm_d)

    def check_imp_type(self):
        """
        Check if something is odd with the force-field parameters (currently
        only eV-units are supported).
        """
        warning_message = (
            "***Warning: Odd improper-type-parameter found. Check if parameter is o.k."
        )

        if self.prm_k > 1 or self.prm_k < 0:
            print(warning_message)

        if int(self.prm_d) not in (-1, 1):
            print(warning_message)


class LongRange(object):
    """
    lr_key: int or str; key of long range interactions
    sitnam_i, sitnam_j: atom names of atom i and atom j
    #TODO: energy units, distance units
    """

    def __init__(
        self,
        lr_key=None,
        atm_key_i=None,
        atm_key_j=None,
        sigma_ij=None,
        epsilon_ij=None,
        acoef=None,
        bcoef=None,
        pairs=None,
        comment=None,
    ):
        """
        pairs:  "ii" or "ij"
        """
        if lr_key is not None:
            self.lr_key = lr_key
        if atm_key_i is not None:
            self.atm_key_i = atm_key_i
        if atm_key_j is not None:
            self.atm_key_j = atm_key_j
        if sigma_ij is not None:
            self.sigma_ij = sigma_ij
        if epsilon_ij is not None:
            self.epsilon_ij = epsilon_ij
        if acoef is not None:
            self.acoef = acoef
        if bcoef is not None:
            self.bcoef = bcoef
        if pairs is not None:
            self.pairs = pairs
        if comment is not None:
            self.comment = comment

    def sig_eps_from_AB(self):
        """
        AB form: A = 4 * eps * sigma**12,
                 B = 4 * eps * sigma**6
                 sigma = (A/B)**(1/6)
                 eps   = B**2/(4*A)
        Sources:    https://en.wikipedia.org/wiki/Lennard-Jones_potential#AB_form
        """
        try:
            self.sigma_ij = (self.acoef / self.bcoef) ** (1 / 6)
        except (ZeroDivisionError):
            self.sigma_ij = 0.00000000e00

        try:
            self.epsilon_ij = self.bcoef ** 2 / (4 * self.acoef)
        except (ZeroDivisionError):
            self.epsilon_ij = 0.00000000e00

        return (self.sigma_ij, self.epsilon_ij)
