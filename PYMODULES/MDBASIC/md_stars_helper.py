from __future__ import print_function
import math
import scipy.constants as sc


# some conversion factors between kj(/mol), kcal(/mol), eV
kCalMol_eV = (1000*sc.calorie)*(1/sc.eV)/sc.N_A
kJ_kCal = 1/sc.calorie
kJMol_eV = 1000*(1/sc.eV)/sc.N_A


def convert_energy_unit(parameter, energy_unit_in, energy_unit_out):
    """
    Unit conversion of epsilon.
    energy_unit_out:   str; 'kcal/mol', 'kj/mol', 'eV'
    """
    if energy_unit_in == "kcal/mol":
        if energy_unit_out == "kj/mol":
            parameter *= sc.calorie
        elif energy_unit_out == "eV":
            parameter *= kCalMol_eV
        else:
            raise KeyError("Wrong keyword for energy_unit_out!")

    elif energy_unit_in == "kj/mol":
        if energy_unit_out == "kcal/mol":
            parameter *= kJ_kCal
        elif energy_unit_out == "eV":
            parameter *= kJMol_eV
        else:
            raise KeyError("Wrong keyword for energy_unit_out!")

    elif energy_unit_in == "eV":
        print("***Warning: eV conversion not implemented (yet)!")
    else:
        raise KeyError("Wrong keyword for energy_unit_out!")

    return parameter


def convert_angle_unit(angle, ang_type_out):
    """
    angle:      float|int; angle to convert
    ang_type:   str; 'rad'|'deg'
    """
    if ang_type_out == "rad":
        return math.radians(angle)
    elif ang_type_out == "deg":
        return math.degrees(angle)
    else:
        print("Wrong ang_type, skipping!")
