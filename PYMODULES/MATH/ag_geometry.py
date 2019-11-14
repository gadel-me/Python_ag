
import math
import itertools as it
import numpy as np
import rmsd
import ag_vectalg as agv
import Transformations as cgt
import pdb

__version__ = "2017-03-30"


#* CENTER OF MASS, DRAW MOLECULAR SYSTEM ***************************************
def get_coord_sys(vt_x, vt_y):
    """
    Return vector z and vector y', so an orthogonal coordinate system
    can be formed.

    Given two vectors x and y, vector b is projected orthogonally
    to vector x (= y').
    The projection (y') and vector x define a plane to which a third
    orthogonal vector z is defined via the cross product. It
    crosses x and y' in the same place.
    From https://en.wikipedia.org/wiki/Scalar_projection
         https://de.wikipedia.org/wiki/Kreuzprodukt
    """
    # The vector projection of a vector 'a' on a nonzero vector 'b' is the
    # orthogonal projection of 'a' onto a straight line parallel to 'b'.
    vt_x = vt_x[:3]  # only xyz-coordinates matter
    vt_y = vt_y[:3]
    # Unit vector ("vt_e") of vector x
    vt_e_x = agv.get_unit_vt(vt_x)
    # Projection of x on y
    vt_y_ex_dot = agv.dot(vt_y, vt_e_x)
    vt_y_ex_dot = [vt_y_ex_dot*i for i in vt_e_x]
    vt_y_proj = [(j - k) for j, k in it.izip(vt_y, vt_y_ex_dot)]
    # Unit vector of vector y
    vt_e_y = agv.get_unit_vt(vt_y_proj)
    # Third orthogonal vector to x and y; vector z (unit vector)
    vt_e_z_ = agv.cross(vt_e_x, vt_e_y)
    vt_e_z = agv.get_unit_vt(vt_e_z_)

    return (vt_e_x, vt_e_y, vt_e_z)


def test_coord_sys(vt_ex, vt_ey, vt_ez, precision=10,
                   coordsys_errorfile="test_coord_sys.err"):
    """
    Test if axis of molecular system made by function 'coord_sys' are orthogonal
    and if each vector is a unit vector.
    """
    coordsys_errorfile = open(coordsys_errorfile, "a")
    # Test vector lengths
    len_ex = round(agv.get_mag(vt_ex), precision)
    len_ey = round(agv.get_mag(vt_ey), precision)
    len_ez = round(agv.get_mag(vt_ez), precision)
    length_warning = "***Warning: {:2.4e} is not a unit vector!\n"

    if len_ex != 1.0:
        coordsys_errorfile.write(length_warning.format(len_ex))
    if len_ey != 1.0:
        coordsys_errorfile.write(length_warning.format(len_ey))
    if len_ez != 1.0:
        coordsys_errorfile.write(length_warning.format(len_ez))

    # Test angles between all three vectors
    exey = round(agv.get_ang(vt_ex, vt_ey), precision)
    exez = round(agv.get_ang(vt_ex, vt_ez), precision)
    eyez = round(agv.get_ang(vt_ey, vt_ez), precision)
    ortho = round(math.pi/2, precision)

    if exey != ortho:
        coordsys_errorfile.write(
            "***Warning: ex and ey are not orthogonal: {}".format(exey)
        )
    if exez != ortho:
        coordsys_errorfile.write(
            "***Warning: ex and ez are not orthogonal: {}".format(exez)
        )
    if eyez != ortho:
        coordsys_errorfile.write(
            "***Warning: ey and ez are not orthogonal: {}".format(eyez)
        )


def get_com(atomic_coordinates, atomic_masses):
    """
    Calculate the center of mass of 1 to n atoms.

    Atomic Masses (u)
    H       1.0079  (DL_POLY: 1.008000)
    C       12.0107 (DL_POLY: 12.011000)
    N       14.0067 (DL_POLY: 14.006700)
    O       15.9994 (DL_POLY: 15.999400)

    From    https://en.wikipedia.org/wiki/Center_of_mass

    Input:
        atomic_coordinates  list;
                            contains all x-, y- and z-atom-coordinates  of
                            which the molecule consists of;
                            the coordinates must be in (1, 3)-lists
        atomic_masses       list;
                            contains all the masses of the atoms in the same
                            order as the coordinates are in 'atomic_coordinates';
    Returns:
        center_of_mass      list;
                            coordinates of center of mass
    """
    # inverse sum of all masses
    # dividing 1 by the sum of masses right away is better for performance
    inv_sum_masses = 1/sum(atomic_masses)
    # multiply each x-, y- and z-coordinate with its corresponding mass
    weighted_coords = []

    for cur_atomic_coord, cur_atomic_mass in it.izip(atomic_coordinates, atomic_masses):
        weighted_coords.append(cur_atomic_coord[0]*cur_atomic_mass)
        weighted_coords.append(cur_atomic_coord[1]*cur_atomic_mass)
        weighted_coords.append(cur_atomic_coord[2]*cur_atomic_mass)

    # sum all x-, y- and z-coordinates up (axis=0 means sum arrays vertically)
    sum_weight_vt = [
        sum(weighted_coords[0::3]),  # x
        sum(weighted_coords[1::3]),  # y
        sum(weighted_coords[2::3])   # z
    ]
    # divide summation of coordinates by the sum of all masses
    center_of_mass = np.array([inv_sum_masses*i for i in sum_weight_vt])
    return center_of_mass


def get_cog(atomic_coordinates):
    """
    Calculate the center of geometry of 1 to n atoms.
    The centre of a polygon is also known as its centroid. It the arithmetic
    mean position of all the points that make up the polygon.

    Sources:    https://en.wikipedia.org/wiki/Centre_(geometry)
                https://deparkes.co.uk/2015/02/28/how-to-find-the-centre-of-a-polygon-in-python/
    """
    num_coords = len(atomic_coordinates)
    center_of_geometry = sum(atomic_coordinates)/num_coords
    return center_of_geometry


#* SPHERE STUFF ****************************************************************
def points_on_sphere(npoints, ndim=3, radius=None):
    """
    From:   http://stackoverflow.com/questions/33976911/generate-a-random-\
            sample-of-points-distributed-on-the-surface-of-a-unit-sphere
            http://mathworld.wolfram.com/SpherePointPicking.html
            https://en.wikipedia.org/wiki/Normal_distribution#\
            Standard_normal_distribution
            M. Muller (1959), doi>10.1145/377939.377946
            https://docs.scipy.org/doc/numpy/reference/generated/\
            numpy.ndarray.transpose.html

    Sample a list of points that are arbitrarily distributed on a sphere.
    Generate a vector consisting of independent samples from three standard
    normal distributions and normalize the vector such its magnitude is 1.
    Vector coordinates are the first numbers of the arrays.

    Input:
        npoints         int;
                        number of points on the sphere that will be generated
        ndim            int;
                        dimension of the vectors (= coordinates of the vectors);
                        2   -> uniformly distributed points on the unit circle
                        3   ->                -"-              the surface of
                        a sphere
                        > 3 ->                -"-              higher dimensions
        radius          float;
                        resize the sphere to 'radius' (= resize length of
                        all vectors to 'radius')
    Returns:
        vec_coords      array;
                        array with 'ndim' sub-arrays (x1, x2, ..., xndim)
                        for the number of dimensions of the coordinates of the
                        vector
    """
    # get coords of the vectors using ndmi=3 for 3-dimensional vectors
    vec_coords = np.random.randn(ndim, npoints)
    # resize all vectors to unit vectors
    vec_coords /= np.linalg.norm(vec_coords, axis=0)

    if radius:
        # resize all vectors to length of radius
        vec_coords *= radius

    # transpose columns to rows
    transposed_coords = vec_coords.transpose()

    return transposed_coords


def point_in_sphere(pt_coords, s_radius, s_origin=[0, 0, 0]):
    """
    FROM:   http://stackoverflow.com/questions/26818772/\
            python-checking-if-a-point-is-in-sphere-with-center-x-y-z
    Check if a point with 'point_coords' is inside a sphere with radius 's_radius'
    and origin 's_origin'.
    """
    in_sphere = False
    px, py, pz = pt_coords
    sx, sy, sz = s_origin

    if (px-sx)**2 + (py-sy)**2 + (pz-sz)**2 < s_radius**2:
        in_sphere = True

    return in_sphere


def point_in_cylinder(pt_coords, c_radius, c_origin=[0, 0, 0]):
    """
    Check if a point (pt_coords) is inside a cylinder with radius 'c_radius'
    with origin 'c_origin'
    """
    in_cylinder = False
    px, py = pt_coords[:2]
    cx, cy = c_origin[:2]

    if (px-cx)**2 + (py-cy)**2 < c_radius**2:
        in_cylinder = True

    return in_cylinder


def arb_rot_matrix(vector):
    """
    Rotate a vector so it aligns with a randomly generated vector. Return the
    rotation matrix.
    """
    vt_u = agv.get_unit_vt(vector)  # scale vector to length of 1
    # generate random vector of length 1
    rand_vt = points_on_sphere(1, ndim=3, radius=1)[0]
    # get angle phi between vt_u and rand_vt
    phi = agv.get_ang(vt_u, rand_vt)
    rot_axis = np.cross(vt_u, rand_vt)
    Mr = cgt.rotation_matrix(phi, rot_axis)
    return Mr


def get_dihedral(ptI, ptJ, ptK, ptL, return_cross=False):
    """
    Define a dihedral/improper between the points I, J, K and L. Defines the planes
    IJK and JKL (if improper is intended to be used, consider right order of atoms).
    Sources: http://cbio.bmt.tue.nl/pumma/index.php/Theory/Potentials
             http://www.vitutor.com/geometry/distance/line_plane.html
             http://kitchingroup.cheme.cmu.edu/blog/2015/01/18/Equation-of-a-plane-through-three-points/
    """
    if not isinstance(ptI, np.ndarray):
        ptI = np.array(ptI)

    if not isinstance(ptJ, np.ndarray):
        ptJ = np.array(ptJ)

    if not isinstance(ptK, np.ndarray):
        ptK = np.array(ptK)

    if not isinstance(ptL, np.ndarray):
        ptL = np.array(ptL)

    # plane IJK
    v1 = ptI - ptJ
    v2 = ptJ - ptK
    # the cross product v1xv2 is a vector normal to both vectors (i.e. to the plane)
    cp1 = np.cross(v1, v2)

    # plane JKL
    v1 = ptJ - ptK
    v2 = ptJ - ptL
    cp2 = np.cross(v1, v2)

    # angle between two planes is the same as the two vectors which are orthogonal
    # to each plane
    if return_cross is True:
        return (cgt.angle_between_vectors(cp1, cp2), cp1, cp2)
    else:
        return cgt.angle_between_vectors(cp1, cp2)


def new_dihedral(p):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.arctan2(y, x)


def get_angle(ptI, ptJ, ptK):
    """
    Bla.

    Get angle between three points I, J and K.
    """
    if not isinstance(ptI, np.ndarray):
        ptI = np.array(ptI)

    if not isinstance(ptJ, np.ndarray):
        ptJ = np.array(ptJ)

    if not isinstance(ptK, np.ndarray):
        ptK = np.array(ptK)

    v1 = ptI - ptJ
    v2 = ptK - ptJ
    return cgt.angle_between_vectors(v1, v2)


def dist_plane_point(plane_vt_1, plane_vt_2, pt):
    """
    Calculate the smallest distance between a point pt and a plane defined
    by vectors plane_vt_1 and plane_vt_2, i.e. the normal between the plane
    and the point.

    Input:
        > plane_vt_1    numpy-array; vector that defines the plane
        > plane_vt_2    numpy-array;            -"-
        > pt            numpy-array; point near the plane

    Returns:
        > vt_d          numpy-array; smallest distance between plane and pt
        > mag           float; length of vt_d
    """
    # get normal vector to plane (vt1 x vt2)
    normal_vt = np.cross(plane_vt_1, plane_vt_2)

    # convert normal vector to a unit vector
    u_normal_vt = np.array(agv.get_unit_vt(normal_vt))

    # get the magnitude between point pt and the plane (= projection of pt on
    # the normal vector)
    mag = np.dot(pt, u_normal_vt)

    # get the vector normal to the plane through point with shortest magnitude
    vt_d = u_normal_vt*mag

    return (vt_d, abs(mag))


def get_molecule_radius(coords, buffering=0):
    """
    Find the largest distance (radius) of each atom towards the center of
    geometry. A buffer may be added to that radius.

    Input:
        > coords    list of numpy-arrays; coordinates of each atom of the molecule
        > buffer    float; value that will be added to the radius in order to enlarge
                    the sphere

    Returns:
        > radius    float; radius of a sphere that envelops the whole molecule
    """

    # center of geometry
    cog = get_cog(coords)
    # initial small radius
    radius = 1e-20

    for coord in coords:
        # current distance
        cdist = np.linalg.norm(cog - coord)
        # replace current radius if current distance is larger
        if cdist > radius:
            radius = cdist

    radius += buffering
    return radius


def get_rmsd(coords, reference):
    """
    Calculate the rmsd of the given coordinates compared to a reference.

    Input:
        > coords        list[np.array, np.array, ...]; coordinates to be compared
                        vs. the reference coordinates
        > reference     list[np.array, np.array, ...]; reference to compare
                        coordinates to

    Sources:   https://pypi.org/project/rmsd/
    """
    # translate center of geometry to the origin for reference and
    # current frame
    coords -= rmsd.centroid(coords)
    reference -= rmsd.centroid(reference)
    return rmsd.kabsch_rmsd(coords, reference)
