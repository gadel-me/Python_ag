from __future__ import print_function, division
import math
import itertools as it

__version__ = "2017-03-30"


#* LINEAR ALGEBRA SECTION ******************************************************
def dot(v1, v2):
    """
    Source: https://en.wikipedia.org/wiki/Dot_product
    Dot-Product in pure python
    """
    return sum(x*y for x, y in it.izip(v1, v2))


def cross(v1, v2):
    """
    Cross product.
    https://en.wikipedia.org/wiki/Cross_product
    """
    a1, a2, a3 = v1
    b1, b2, b3 = v2
    return (a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1)


def mm_mult(m1, m2):
    """
    Source: http://www.programiz.com/python-programming/examples/multiply-matrix
    Multiplication of two matrices of same size
    Input:
        m1      array of arrays
        m2      array of arrays
    """
    result = [[sum(a*b for a, b in it.izip(X_row, Y_col)) for Y_col in it.izip(*m2)] for
              X_row in m1]
    return result


def mv_mult(m, v):
    """
    Multiplication of matrix with vector.
    """
    if isinstance(v, tuple):
        v = list(v)

    if len(v) == 3:
        v.append(1)

    return [dot(row, v) for row in m]


def vs_mult(vector, skalar):
    """
    Multiply a vector with a skalar.
    """
    vt_multiplied = [i*skalar for i in vector]
    return vt_multiplied


def get_vt(tail, head, direction="head"):
    """
    Return vector defined by 2 Points (tail(x/y/z/.../n, head(x'/y'/z'/.../n)).
    Source: http://mathworld.wolfram.com/Vector.html

    Input:
        tail:       list;
                    may contain n-elements
        head:       list;
                    may contain n-elements
        direction:  str;
                    'head' or 'tail' allowed;
                    defines direction towards v_out will point
    Returns:
        list

    """
    v_out = [j-i for i, j in it.izip(tail, head)]

    if direction == "head":
        pass
    elif direction == "tail":
        v_out = [-1*i for i in v_out]
    else:
        raise RuntimeError("{} is not a proper direction".format(direction))

    return v_out


def add_vts(*vects):
    """
    Add all given vectors.
    Only works with vector coordinates as lists.
    """
    return [sum(x) for x in it.izip(*vects)]


def get_mag(vector, unrooted=False):
    """
    Check magnitude of a vector fast, i.e. not using sqrt (root is slow!)
    Source: http://mathworld.wolfram.com/Norm.html
    """
    _premagn = [i**2 for i in vector]
    _fastmagn = sum(_premagn)

    if not unrooted:
        return math.sqrt(_fastmagn)
    else:
        return _fastmagn


def get_unit_vt(vector):
    """
    A unit vector is a vector of length 1, sometimes also called
    a direction vector (Jeffreys and Jeffreys 1988).
    Convert a given vector to its unit vector.
    Source: http://mathworld.wolfram.com/UnitVector.html
    """
    _magn = get_mag(vector, unrooted=False)
    unit_vector = [(1/_magn)*i for i in vector]
    return unit_vector


def get_ang(v1, v2, deg=False):
    """
    Input:
        v1        list;
                  vector 1
        v2        list;
                  vector 2
    Returns the angle between vectors 'v1' and 'v2':

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966 (1/2 pi)
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793 (pi)
    """
    numerator = sum([i*j for i, j in it.izip(v1, v2)])
    denominator = get_mag(v1)*get_mag(v2)
    gamma = math.acos(numerator/denominator)

    # return angle in degrees instead of radians
    if deg:
        gamma = math.degrees(gamma)

    return gamma


def get_plane(vt_1, vt_2, vt_p=[0, 0, 0]):
    """
    Define a plane by the vectors vt_1, vt_2 and the positional vector vt_p.
    Calculate the normal vector of vt_1 and vt_2 (which defines a, b, c) and
    the dot product of vt_p with vt_n which defines d.

    Returns the variables a, b, c, d for the plane equation: ax + by + cz +d = 0
    which checks if a point (x, y, z) lies on the plane.

    Sources: https://en.wikipedia.org/wiki/Plane_(mathematics)
             https://de.wikipedia.org/wiki/Koordinatenform

    Input:
        > vt_1   list; first vector to define the plane
        > vt_2   list; second vector to define the plane
        > vt_p  list; position-vector p

    Returns:
        a, b, c, d      float;
    """
    vt_n = get_unit_vt(cross(vt_1, vt_2))  # normal vector (perpendicular to vt_1 and vt_2)
    a, b, c = vt_n
    d = dot(vt_p, vt_n)
    return (a, b, c, d)


def get_point_plane_dist(p, a, b, c, d, distance_sign_only=False):
    """
    Calculate the distance from a point p to the plane defined by a, b, c and d.
    Sources: https://en.wikipedia.org/wiki/Plane_(mathematics)#Describing_a_plane_through_three_points

    Input:
        > p                    list; coordinates of the point
        > a, b, c, d           float; plane variables
        > distance_sign_only   boolean; just solve the plane equation and return
                               a positive or negative float (just to check on which
                               side of the cube the point is)
    """
    if distance_sign_only is True:
        distance = a*p[0]+b*p[1]+c*p[2]-d
    else:
        distance = abs(a*p[0]+b*p[1]+c*p[2]-d)/math.sqrt(a**2+b**2+c**2)

    return distance
