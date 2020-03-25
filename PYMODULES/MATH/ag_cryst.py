import math
import numpy as np
#TODO from typing import Union, Tuple, List
import ag_vectalg as agv

__version__ = "2017-03-30"


# * CELLULAR BOX SECTION; BOX MEANS (UNIT/SUPER-)CELL
def box_lat_volume(a, b, c, alpha, beta, gamma):
    """
    Calculate the volume of a cell which is defined by its vectors and angles.
    From: https://en.wikipedia.org/wiki/Fractional_coordinates
    Input:
        a           float;
                    cell vector a
        b           float;
                    cell vector b
        c           float;
                    cell vector c
        alpha       float;
                    cell angle alpha in radians
        beta        float;
                    cell angle beta in radians
        gamma       float;
                    cell angle gamma in radians
    """
    v_1 = 1 - math.cos(alpha) ** 2 - math.cos(beta) ** 2 - math.cos(gamma) ** 2
    v_2 = 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma)
    v_3 = a * b * c
    volume = v_3 * math.sqrt(v_1 + v_2)
    return volume


def box_cart2lat(vector_a, vector_b, vector_c):
    """
    Convert box vectors
    Cartesian -> Lattice
    """
    a = agv.get_mag(vector_a)
    b = agv.get_mag(vector_b)
    c = agv.get_mag(vector_c)
    alpha = agv.get_ang(vector_b, vector_c)
    beta = agv.get_ang(vector_a, vector_c)
    gamma = agv.get_ang(vector_a, vector_b)
    return [a, b, c, alpha, beta, gamma]


def box_lat2cart(a, b, c, alpha, beta, gamma):
    """
    Convert box vectors
    Lattice -> Cartesian
    """
    M_f2c = M_fract2cart(a, b, c, alpha, beta, gamma)
    # to get right handed box
    a = [1, 0, 0]
    b = [0, 1, 0]
    c = [0, 0, 1]
    vector_a = agv.mv_mult(M_f2c, a)
    vector_b = agv.mv_mult(M_f2c, b)
    vector_c = agv.mv_mult(M_f2c, c)
    return (vector_a, vector_b, vector_c)


def box_lat2lmp(a, b, c, alpha, beta, gamma):
    """
    Lattice -> Lammps
    From:   http://lammps.sandia.gov/doc/Section_howto.html#howto-12
    """
    # (lx,ly,lz) = (xhi-xlo,yhi-ylo,zhi-zlo) and tilt factors (xy,xz,yz) is as follows:
    lx = a
    xy = b * math.cos(gamma)
    xz = c * math.cos(beta)
    ly2 = b ** 2 - xy ** 2
    ly = math.sqrt(ly2)
    yz = (b * c * math.cos(alpha) - xy * xz) / (ly)
    lz2 = c ** 2 - xz ** 2 - yz ** 2
    lz = math.sqrt(lz2)

    return (lx, ly, lz, xy, xz, yz)


def box_cart2lmp(vector_a, vector_b, vector_c):
    """
    Cartesian -> Lammps
    From:   http://lammps.sandia.gov/doc/Section_howto.html#howto-12
    """
    lx = vector_a[0]
    ly = vector_b[1]
    lz = vector_c[2]

    xy = vector_b[0]
    xz = vector_c[0]
    yz = vector_c[1]

    return (lx, ly, lz, xy, xz, yz)


def box_lmp2lat(lx, ly, lz, xy, xz, yz):
    """
    Lammps -> Lattice
    From:   http://lammps.sandia.gov/doc/Section_howto.html#howto-12
    """
    a = lx
    b = math.sqrt(ly ** 2 + xy ** 2)
    c = math.sqrt(lz ** 2 + xz ** 2 + yz ** 2)

    cos_alpha = (xy * xz + ly * yz) / (b * c)
    cos_beta = xz / c

    cos_gamma = xy / b
    alpha = math.acos(cos_alpha)
    beta = math.acos(cos_beta)
    gamma = math.acos(cos_gamma)

    return (a, b, c, alpha, beta, gamma)


def box_lmp2cart(lx, ly, lz, xy, xz, yz):
    """
    Lammps -> Cartesian
    Gives right handed coordinate system.
    """
    vector_a = [lx, 0, 0]
    vector_b = [xy, ly, 0]
    vector_c = [xz, yz, lz]

    return (vector_a, vector_b, vector_c)

def box_lmp2alat(lx, ly, lz, xy, xz, yz):
    alat = lx
    vector_a = [lx / alat, 0, 0]
    vector_b = [xy / alat, ly / alat, 0]
    vector_c = [xz / alat, yz / alat, lz / alat]

    return (vector_a, vector_b, vector_c, alat)

def box_cart2alat(vector_a, vector_b, vector_c):
    """box_cart2alat [summary]

    Parameters
    ----------
    vector_a : Union[list, tuple]
        [description]
    vector_b : Union[list, tuple]
        [description]
    vector_c : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    # check if vector a is aligned to x-axis
    alignement_err = "Vector a is not aligned to the x-axis!"
    assert vector_a[1] == 0, alignement_err
    assert vector_a[2] == 0, alignement_err

    # TODO: Rotate box or return rotational matrix at least so the alignment
    # TODO: is fine (vector_a => x-axis)

    # calculate magnitude of vector_a
    alat = agv.get_mag(vector_a)
    # convert vector_a to unit vector
    vector_a = [i / alat for i in vector_a]
    vector_b = [i / alat for i in vector_b]
    vector_c = [i / alat for i in vector_c]
    return (vector_a, vector_b, vector_c, alat)

# * FRACTIONAL COORDINATES SECTION **********************************************
def M_fract2cart(a, b, c, alpha, beta, gamma):
    """
    From:   https://en.wikipedia.org/wiki/Fractional_coordinates
    """
    a11 = a
    a12 = 0
    a13 = 0

    a21 = b * math.cos(gamma)
    a22 = b * math.sin(gamma)
    a23 = 0

    a31 = c * math.cos(beta)
    a32_enum = math.cos(alpha) - (math.cos(beta) * math.cos(gamma))
    a32_denom = math.sin(gamma)
    a32 = c * (a32_enum / a32_denom)

    a33_enum = box_lat_volume(a, b, c, alpha, beta, gamma)
    a33_denom = a * b * math.sin(gamma)
    a_33 = a33_enum / a33_denom

    return [[a11, a21, a31], [a12, a22, a32], [a13, a23, a_33]]


def M_cart2fract(a, b, c, alpha, beta, gamma):
    """
    From:   https://en.wikipedia.org/wiki/Fractional_coordinates
    """
    vol = box_lat_volume(a, b, c, alpha, beta, gamma)

    # column 1
    a11 = 1 / a
    a12 = 0
    a13 = 0

    # column 2
    a21_enum = math.cos(gamma)
    a21_denom = a * math.sin(gamma)
    a21 = -1 * (a21_enum / a21_denom)

    a22_enum = 1
    a22_denom = b * math.sin(gamma)
    a22 = a22_enum / a22_denom

    a23 = 0

    # column 3
    a31_enum = math.cos(alpha) * math.cos(gamma) - math.cos(beta)
    a31_denom = vol * math.sin(gamma)
    a31 = b * c * (a31_enum / a31_denom)

    a32_enum = math.cos(beta) * math.cos(gamma) - math.cos(alpha)
    a32_denom = vol * math.sin(gamma)
    a32 = a * c * (a32_enum / a32_denom)

    a33_enum = math.sin(gamma)
    a33_denom = vol
    a_33 = a * b * (a33_enum / a33_denom)

    return [[a11, a21, a31], [a12, a22, a32], [a13, a23, a_33]]


# debugging area
if __name__ == '__main__':
    a,b,c, alat = box_cart2alat([4, 1, 0], [2, 2, 2], [4, 4, 6])
