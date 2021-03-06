#!/usr/bin/env python

import pdb
import os
import time
import numpy as np
from collections import OrderedDict
from pyquaternion import Quaternion
from PIL import Image
import ag_unify_md as agum
import ag_geometry
import ag_vectalg as agv
import Transformations as cgt
import pantone_colors

# VMD MODULES
import atomsel   # Replaces AtomSel and atomselection
import axes
import color
import display
import graphics
import imd
import label
import material
import molecule
import molrep
import mouse
import render
import trans
import vmdmenu
import Label
import Material
import Molecule
import VMD


#=============================================================================================#
# VMD HELPER FUNCTIONS
#=============================================================================================#
PANTONE_COLORS = pantone_colors.pantone_colors_hex


def vmd_draw_angle(molid, frame, atomids=None, atm_coords=None, canvas=False, resolution=200, radius=0.5, drawcolor="blue", cylinder_radius=0.01):
    """
    Draws part of a circle to visualize the measured angle.

        Input:
            > moldid     int; index of molecule to draw the angle for
            > frame      int; frame number to draw the angle for
            > atomids    tuple; all three atom ids with second one as angular
                         point
            > canvas     boolean; draw a canvas between bow and angle axis
            > resolution int; number of cylinders and spheres that form the angle

    """
    if atomids is not None:
        id1, id2, id3 = atomids
        # select all atoms
        sel = atomsel.atomsel("all", molid, frame)
        # get coordinates
        x = sel.get("x")
        y = sel.get("y")
        z = sel.get("z")
        # get rotational axis
        atm1Crds = np.array([x[id1], y[id1], z[id1]])
        atm2Crds = np.array([x[id2], y[id2], z[id2]])
        atm3Crds = np.array([x[id3], y[id3], z[id3]])
    elif atm_coords is not None:
        atm1Crds, atm2Crds, atm3Crds = atm_coords
    else:
        raise Warning("At least atomic coordinates or atom ids have to be provided!")

    # draw bows in the current molecule, transparent canvas in a new molecule
    if canvas is True:
        molecule.new("angle canvas")
        graphics_id = molecule.num() - 1  # ids start with 0
    else:
        graphics_id = molid

    graphics.color(graphics_id, drawcolor)

    #pdb.set_trace()
    # define angle and rotational axis
    vt1    = atm1Crds - atm2Crds
    vt2    = atm3Crds - atm2Crds
    vtRot  = np.cross(vt1, vt2)

    # unit vectors for quaternions
    vt1   /= np.linalg.norm(vt1)
    vt2   /= np.linalg.norm(vt2)
    vtRot /= np.linalg.norm(vtRot)

    # define angle offset
    max_angle = np.arccos(np.dot(vt1, vt2))  # denominator is 1 since vectors are normed
    print("Measured angles: {}".format(np.degrees(max_angle)))
    angle_offset = np.pi / resolution

    vt1 *= radius
    vt2 *= radius

    for cntr, theta in enumerate(np.arange(0, max_angle, angle_offset)):
        # previous vector
        if cntr == 0:
            vt_pre = vt1 + atm2Crds

        # define quaternion
        q = Quaternion(axis=vtRot, angle=theta)
        #print(q)

        # rotate vector
        vt_cur = q.rotate(vt1)
        vt_cur += atm2Crds

        if canvas is True and cntr != 0:
            graphics.triangle(graphics_id, tuple(vt_pre), tuple(atm2Crds), tuple(vt_cur))
        else:
            graphics.cylinder(graphics_id, tuple(vt_pre), tuple(vt_cur), radius=cylinder_radius,
                              resolution=20, filled=1)
        vt_pre = vt_cur

    if canvas is True:
        graphics.triangle(graphics_id, tuple(vt_pre), tuple(atm2Crds), tuple(vt2 + atm2Crds))
        graphics.material(graphics_id, "Transparent")
    else:
        graphics.cylinder(graphics_id, tuple(vt_pre), tuple(vt2 + atm2Crds), radius=cylinder_radius,
                          resolution=20, filled=1)


def draw_torsion(molid, frame, atomids, canvas=False, resolution=50, radius=0.5, drawcolor="blue"):
    """
    """
    id1, id2, id3, id4 = atomids

    # draw bows in the current molecule, transparent canvas in a new molecule
    if canvas is True:
        molecule.new("angle canvas")
        graphics_id = molecule.num() - 1  # ids start with 0
    else:
        graphics_id = molid

    graphics.color(graphics_id, drawcolor)

    # select all atoms
    sel = atomsel.atomsel("all", molid, frame)

    # get coordinates
    x = sel.get("x")
    y = sel.get("y")
    z = sel.get("z")

    # get rotational axis
    atm1Crds = np.array([x[id1], y[id1], z[id1]])
    atm2Crds = np.array([x[id2], y[id2], z[id2]])
    atm3Crds = np.array([x[id3], y[id3], z[id3]])
    atm4Crds = np.array([x[id4], y[id4], z[id4]])


    # define angle and rotational axis
    #vt1    = atm1Crds - atm2Crds
    #vt2    = atm3Crds - atm2Crds
    #vtRot  = np.cross(vt1, vt2)

    # unit vectors for quaternions
    #vt1   /= np.linalg.norm(vt1)
    #vt2   /= np.linalg.norm(vt2)
    #vtRot /= np.linalg.norm(vtRot)
    #print(atm1Crds)
    _, cp1, cp2 = ag_geometry.get_dihedral(atm1Crds, atm2Crds, atm3Crds, atm4Crds, return_cross=True)

    # draw vectors between id2 and id3, therefor get the relevant vector
    center_pt = atm2Crds + (atm3Crds - atm2Crds) * 0.5
    cp1 += center_pt
    cp2 += center_pt
    vmd_draw_angle(molid, frame, atm_coords=[cp1, center_pt, cp2], canvas=False, resolution=resolution, radius=radius, drawcolor="blue")


def vmd_label(molid, key, selection="all", label_color="white", textsize=1.0, offset=(0.0, 0.0, 0.0), idx_offset=None):
    """
    Labels atoms by index, charge, etc.

        Input:
           > molid       int; index of molecule to label
           > key         str; attribute to label atom by
           > textsize    float; size of label font
           > offset      (1,3)-list; offset to shift the label by
           > selection   str; selected atoms
           > num_offset  int; offset to add to atomic index
    """

    # select atoms to label
    selected_atoms = atomsel.atomsel(selection, molid)

    # complete selection of all atoms in order to label just the ones needed
    complete_selection = atomsel.atomsel("all", molid)
    atom_idxs = selected_atoms.get("index")
    print(len(atom_idxs))

    # get its coordinates
    xs = complete_selection.get("x")
    ys = complete_selection.get("y")
    zs = complete_selection.get("z")
    values = complete_selection.get(key)

    if key == "index":
        if idx_offset is not None:
            values = [i + idx_offset for i in values]

    if key == "charge":
        total_charge = sum(values)
        print("Total charge is {}".format(total_charge))
        del total_charge

    # set a label color
    graphics.color(molid, label_color)

    for atom_idx in atom_idxs:
        label_pos = (xs[atom_idx] + offset[0], ys[atom_idx] + offset[1], zs[atom_idx] + offset[2])

        if key == "index":
            labelformat = "{}"
        else:
            labelformat = "{:2.3f}"

        label_text = labelformat.format(values[atom_idx])
        graphics.text(molid, tuple(label_pos), label_text, textsize)


def vmd_draw_arrow(molid, start, end, cylinder_radius=0.4, cone_radius=1.0,
                   cone_length=0.15, resolution=50, double_arrow=False,
                   vmd_material="Basic1Pantone", drawcolor="blue", lstyle="solid"):
    """
    Draws an arrow from start to end using the arrow color and a certain radius
    for the cylinder and the cone (top).
    Draw cylinder of 90 % of the max length and the tip from 75 % of the total
    arrow length.

    Input:
        > molid             int; id of the molecule to draw the arrow to
        > start             np-array; (1,3)-tuple with starting coordinates
        > end               np-array; (1,3)-tuple with ending coordinates
        > cylinder_radius   float; radius of the arrow base
        > cone_radius       float; radius of the cone

    lstyle : str; solid | dashed
        draw a solid (default) or dashed arrow
    """
    graphics.material(molid, vmd_material)
    graphics.color(molid, drawcolor)
    p_vector = end - start

    # shift starting point of the vector by cone_length towards the end point
    if double_arrow is True:
        # shorten vector so there is place for the cone at the starting point
        p_start = start + p_vector * cone_length
    else:
        p_start = start

    # shorten vector so there is place for the cone at the ending point
    p_end = start + p_vector * (1 - cone_length)
    graphics.cone(molid, tuple(p_end), tuple(end), radius=cone_radius, resolution=resolution)
    graphics.cone(molid, tuple(p_start), tuple(start), radius=cone_radius, resolution=resolution)

    if lstyle == "solid":
        graphics.cylinder(molid, tuple(p_start), tuple(p_end), radius=cylinder_radius, resolution=resolution)
    elif lstyle == "dashed":
        # length of a single cylinder
        len_cylinder = 0.08

        # normed length of the vector between the cones
        p_vector_len = np.linalg.norm(p_start - p_end)

        # number of cylinders that fit between the cones with len_cylinder length
        ncylinders = int(p_vector_len // len_cylinder)
        p_vector_unit = (p_end - p_start) / p_vector_len

        # define the intermediate vector
        inter_end = [0, 0, 0]

        # draw vectors ncylinders times
        for i in range(ncylinders):
            inter_end = p_start + p_vector_unit * len_cylinder

            if i % 2 == 0:
                graphics.cylinder(molid, tuple(p_start), tuple(inter_end), radius=cylinder_radius, resolution=resolution)
                #pdb.set_trace()

            p_start = inter_end
    else:
        pass


def vmd_render_scene(image_out, image_size=[2000, 2000], renderer="TachyonLOptiXInternal", frame_idx=0):
    """
    Render the image using renderer with resolution image_size.
    """
    disp_default_size = display.get("size")
    # higher resolution
    display.set(size=image_size)
    display.update()
    display.update_ui()
    VMD.evaltcl("light 0 off")
    render.render(renderer, "{}.ppm".format(image_out))
    VMD.evaltcl("light 0 on")
    image = Image.open("{}.ppm".format(image_out))
    image.save("{}.png".format(image_out), format="PNG")
    #display.set(size=disp_default_size)
    display.update()
    display.update_ui()
    os.remove("{}.ppm".format(image_out))


def vmd_prepare_scene():
    """Change the lighting and background so it is the same for each image when rendering."""
    VMD.evaltcl("light 1 off")
    VMD.evaltcl("display shadows off")
    color_display = color.get_colormap("Display")
    color_display["Background"] = "white"
    color.set_colormap("Display", color_display)
    display.set(depthcue=0)
    display.update()
    display.update_ui()


def vmd_load_molecule(coordsfile, filetype="lammpsdata", dcd=None,
                      selection="all", scale=1.0, molid=0,
                      vmd_material="Basic1Pantone", style="CPK 1.000000 0.300000 12.000000 12.000000"):
    """
    """
    if not molecule.exists(molid):

        # load molecule in vmd
        if dcd is not None:
            molecule.load(filetype, coordsfile, "dcd", dcd)
        else:
            molecule.load(filetype, coordsfile)

        molrep.modrep(molid, 0, style=style, material=vmd_material, color="Name", sel=selection)
        #trans.scale(scale)
        VMD.evaltcl("scale to {}".format(scale))


def vmd_draw_box(lines, molid=0, vmd_material="Basic1Pantone", drawcolor="blue",
                 lstyle="solid", lwidth=1, lfreq=15):
    """
    Draw the unit cell in molid.

    This function takes the output from the lines_to_draw function.

    Parameters
    ----------
    molid : int, optional
        id of loaded molecule in vmd

    lines : {list{list, list}, list{list, list}}
        all lines with starting and endpoints that are to be drawn;
        this is the output from the lines_to_draw function

    lstyle : str; 'solid' or 'dashed' or '--'
        draw dashed lines or own dashed lines

    lwidth : int
        width of the line

    lfreq : int
        frequency a line shall be drawn, which is lfreq / 2 since every
        second line is skipped

    """
    #all_molids = molecule.listall()

    graphics.material(molid, vmd_material)
    graphics.color(molid, drawcolor)

    for ucell_line in lines:
        #print(ucell_line)

        # use a more rough dashing than the built-in one
        if lstyle == "--":
            # divide current line into sublines
            pstart = ucell_line[0]

            # get the subline and norm it
            subline = ucell_line[1] - ucell_line[0]
            subline /= lfreq

            for cntr in range(lfreq):
                pstop = pstart + subline

                if cntr % 2 == 0:
                    #print(cntr)
                    graphics.line(molid, tuple(pstart), tuple(pstop),
                                  style="solid", width=lwidth)

                pstart = pstop
        else:
            graphics.line(molid, tuple(ucell_line[0]), tuple(ucell_line[1]),
                          style=lstyle, width=lwidth)


def vmd_draw_lattice_axes(ucell_a, ucell_b, ucell_c, origin=(-4, -4, -4),
                          molid=0, vmd_material="Basic1Pantone", label_axis=True, textsize=1.0,
                          offset=(0.0, 0.0, 0.0)):
    """
    Draw the arrows that show the coordinate system of the box according to the lattice vectors.

    Parameters
    ----------
    ucell_a : np-array
        cell vector a

    ucell_b : np-array
        cell vector b

    ucell_c : np-array
        cell vector c

    origin : tuple or list or numpy.ndarray {float float float}
        origin where to place the lattice axis

    molid : int; default = 0
        molid in vmd to draw the lattice axes in

    vmd_material : str
        the material to draw the lattice axes in

    offset : tuple or list or numpy.ndarray {float float float}
        offset where to place the letters for the lattice coordinates

    """
    if not isinstance(origin, np.ndarray):
        origin = np.array(origin)

    if not isinstance(offset, np.ndarray):
        offset = np.array(offset)

    # change vector size to 8
    ucell_a = ucell_a / np.linalg.norm(ucell_a) * 8
    ucell_b = ucell_b / np.linalg.norm(ucell_b) * 8
    ucell_c = ucell_c / np.linalg.norm(ucell_c) * 8

    graphics.material(molid, vmd_material)

    # axis a
    vmd_draw_arrow(molid, origin, origin + ucell_a, cylinder_radius=0.4, cone_radius=1.0, drawcolor="red")
    label_pos = origin + ucell_a * 1.05 + offset

    if label_axis is True:
        graphics.text(molid, tuple(label_pos), "a", textsize)

    # axis a
    vmd_draw_arrow(molid, origin, origin + ucell_b, cylinder_radius=0.4, cone_radius=1.0, drawcolor="green")
    label_pos = origin + ucell_b * 1.05 + offset

    if label_axis is True:
        graphics.text(molid, tuple(label_pos), "b", textsize)

    # axis c
    vmd_draw_arrow(molid, origin, origin + ucell_c, cylinder_radius=0.4, cone_radius=1.0, drawcolor="blue")
    label_pos = origin + ucell_c * 1.05 + offset

    if label_axis is True:
        graphics.text(molid, tuple(label_pos), "c", textsize)


# FUNCTIONS TO PLOT SUPERCELLS OF CRYSTAL STRUCTURES
def ucell_scell_factors(lmpdat_ucell, lmpdat_scell):
    """
    Get the factors the unit cell was multiplied by to gain the super cell.

    Parameters
    ----------
    lmpdat_ucell : str
        lammps data file of the unit cell

    lmpdat_scell : str
        lammps data file of the super cell

    """
    supercell = agum.Unification()
    supercell.read_lmpdat(lmpdat_scell)
    unitcell = agum.Unification()
    unitcell.read_lmpdat(lmpdat_ucell)
    supercell.ts_boxes[-1].box_lmp2lat()
    unitcell.ts_boxes[-1].box_lmp2lat()
    fa = int(supercell.ts_boxes[-1].ltc_a / unitcell.ts_boxes[-1].ltc_a)
    fb = int(supercell.ts_boxes[-1].ltc_b / unitcell.ts_boxes[-1].ltc_b)
    fc = int(supercell.ts_boxes[-1].ltc_c / unitcell.ts_boxes[-1].ltc_c)
    del (supercell, unitcell)
    return [fa, fb, fc]


def lines_to_draw(p0, ucell_a, ucell_b, ucell_c):
    """
    Get the starting and end points of each line to draw for the unit cell.

    Parameters
    ----------
    p0 : np-array
        origin where to place the line beginning

    ucell_a : np-array
        cell vector a

    ucell_b : np-array
        cell vector b

    ucell_c : np-array
        cell vector c

    """
    vectors = []

    # vector a
    for i in [0, 1]:
        for j in [0, 1]:
            vstart = p0 + i * ucell_b + j * ucell_c
            vstop = vstart + ucell_a
            cvect = [vstart, vstop]
            vectors.append(cvect)

    # vector b
    for i in [0, 1]:
        for j in [0, 1]:
            vstart = p0 + i * ucell_a + j * ucell_c
            vstop = vstart + ucell_b
            cvect = [vstart, vstop]
            vectors.append(cvect)

    # vector c
    for i in [0, 1]:
        for j in [0, 1]:
            vstart = p0 + i * ucell_a + j * ucell_b
            vstop = vstart + ucell_c
            cvect = [vstart, vstop]
            vectors.append(cvect)

    return vectors


def ucell_vects(supercell_dcd, fa, fb, fc, frame_id=-1):
    """
    Get the unit cell vectors from a super cell according to the multiplication factors.

    Parameters
    ----------
    supercell_dcd : str
        the dcd file to get the super cell lattice vectors from

    fa : int or float
        the factor to divide the super cell lattice vector a by in order
        to get the unit cell vector a

    fb : int or float
        the factor to divide the super cell lattice vector b by in order
        to get the unit cell vector b

    fc : int or float
        the factor to divide the super cell lattice vector c by in order
        to get the unit cell vector c

    """
    scell = agum.Unification()
    scell.import_dcd(supercell_dcd)
    scell.read_frames()
    scell.ts_boxes[frame_id].box_lmp2cart()
    ucell_a = np.array(scell.ts_boxes[frame_id].crt_a) / fa
    ucell_b = np.array(scell.ts_boxes[frame_id].crt_b) / fb
    ucell_c = np.array(scell.ts_boxes[frame_id].crt_c) / fc
    return (ucell_a, ucell_b, ucell_c)


def draw_dotted_line(molid, start, end, sradius=0.01, drawcolor="blue", sdist=0.1, vmd_material="Basic1Pantone"):
    """
    Draw a dotted line between two pints 'start' and 'end'.

    Parameters
    ----------
    molid : int
        id of the molecule to draw the arrow to

    start : np-array; (1,3)-tuple
        starting coordinates

    end : np-array; (1,3)-tuple
        ending coordinates

    sradius : float
        radii of the dots

    drawcolor : str
        color of the dots

    sdist : float
        distance between points

    """
    graphics.material(molid, vmd_material)
    graphics.color(molid, drawcolor)

    # vector between start and end
    p_vector = end - start
    length_p_vector = np.linalg.norm(p_vector)

    # normed p_vector
    normed_p_vector = p_vector / length_p_vector

    # number of points
    num_points = int(length_p_vector / sdist)

    # first sphere
    for pt in range(num_points):
        csphere = start + normed_p_vector * sdist * pt
        graphics.sphere(molid, center=tuple(csphere), radius=sradius)


def vmd_coords_numpy_arrays(molid, frame_id, selection="all"):
    """
    """
    sel = atomsel.atomsel(selection, molid, frame_id)
    # get coordinates
    x = sel.get("x")
    y = sel.get("y")
    z = sel.get("z")
    coords = [np.array([i, j, k]) for i, j, k in zip(x, y, z)]

    # omit the list the coords are in, if the coordinates of only one atom are desired
    if len(coords) == 1:
        return coords[0]
    else:
        return coords


def draw_double_bonds(molid, frame_id, bonds, radius=0.02, filled=True, vmd_material="Basic1Pantone"):
    """
    Draw double bonds using cylinders given a dictionary with according bonds.

    Parameters
    ----------
    molid : int
        molecule id to draw the cylinder in
    frame_id : int
        id of the frame to use the coordinates from
    bonds : list {list {int, int, int, int}, list {int, int, int, int}, ...}
        first two integers of each sublist form the bond, the others form
        the vector which forms the plane in which the double will be located
    """

    # TODO: since a double bond is always in a plane, both neighbor bonds
    # TODO: are also in that plane, therefor a plane can always easily
    #       be defined without explicitly naming the neighbors
    selection_str = "index {}"
    element_colors = color.get_colormap("Element")
    graphics.material(molid, vmd_material)

    for cbond in bonds:
        selection_str1 = selection_str.format(cbond[0])
        selection_str2 = selection_str.format(cbond[1])
        selection_str3 = selection_str.format(cbond[2])
        selection_str4 = selection_str.format(cbond[3])

        # get vector which is the bond
        vt1_coords = vmd_coords_numpy_arrays(molid, frame_id, selection_str1)
        vt2_coords = vmd_coords_numpy_arrays(molid, frame_id, selection_str2)
        #print(vt1_coords)
        bond = vt2_coords - vt1_coords

        # get vector forms the plane with bond vector
        vt3_coords = vmd_coords_numpy_arrays(molid, frame_id, selection_str3)
        vt4_coords = vmd_coords_numpy_arrays(molid, frame_id, selection_str4)
        bond_neigh = vt4_coords - vt3_coords

        # get normal to bond vector and neighbor vector
        normal1 = np.cross(bond_neigh, bond)
        normal1 /= np.linalg.norm(normal1)

        # get normal to normal of bond and bond neigh
        # and bond - this is the desired shift vector for the bond
        normal2 = np.cross(bond, normal1)
        normal2 /= np.linalg.norm(normal2)

        # draw double bond
        tail_bond2 = vt1_coords + normal2 * 0.1
        head_bond2 = vt2_coords + normal2 * 0.1

        # center between tail and head (for coloring)
        center = bond / 2 + vt1_coords + normal2 * 0.1

        # get color of atoms
        selection1 = atomsel.atomsel(selection_str1)
        selection2 = atomsel.atomsel(selection_str2)

        # get element name
        atom1_element = selection1.get("element")[0]
        atom2_element = selection2.get("element")[0]

        # get element color
        atom1_color = element_colors[atom1_element]
        atom2_color = element_colors[atom2_element]

        # draw double bond - atom 1 to center using corresponding color
        graphics.color(molid, atom1_color)
        graphics.cylinder(molid, tuple(tail_bond2), tuple(center), radius=radius, filled=filled)

        # draw double bond - center to atom 2 using corresponding color
        graphics.color(molid, atom2_color)
        graphics.cylinder(molid, tuple(center), tuple(head_bond2), radius=radius, filled=filled)


def draw_circle(molid, rotaxis, center, angle=360, drawcolor="blue", circle_radius=0.2,
                vmd_material="Basic1Pantone", cylinder_radius=0.01,
                cone_radius=0.02, cone_length=10):
    """
    Draw a circle around center using a rotational axis.

    Parameters
    ----------
    molid : int
    frame_id : int

    rotaxis : numpy-array {float float float}
        rotational axis

    center : numpy-array {float float float}
        center to rotate around

    angle : float or int
        angle to rotate about in degrees

    color : str
        color of the circle

    circle_radius : float
        circle_radius of the circle

    resolution : int
        vectors the circle consists of

    """
    # change color and material of circle
    graphics.color(molid, drawcolor)
    graphics.material(molid, vmd_material)

    # norm rotational axis
    rotaxis /= np.linalg.norm(rotaxis)

    # create and scale first vector according to given radius
    normal1 = np.cross(rotaxis, center)
    normal1 /= np.linalg.norm(normal1)
    normal1 *= circle_radius

    vector_start = normal1

    # double angle for more points
    angle *= 2

    for cangle in range(angle - cone_length):
        #print(cangle)
        # convert angle to radians
        cangle = np.radians(cangle / 2)
        #print(cangle)

        # define quaternion according to current angle
        quaternion = Quaternion(axis=rotaxis, angle=cangle)

        # define end of vector
        vector_end = quaternion.rotate(normal1)
        #vector_end /= np.linalg.norm(vector_end)
        #vector_end *= radius

        graphics.cylinder(molid,
                          tuple(center + vector_start),
                          tuple(center + vector_end),
                          radius=cylinder_radius, resolution=20, filled=1)
        #graphics.sphere(molid, tuple(center + vector_start), radius=0.02)

        # new starting point is old ending
        vector_start = vector_end

    # draw an arrow beginning 5 degrees before the last angle
    vector_start = vector_end

    quaternion = Quaternion(axis=rotaxis, angle=np.radians(angle / 2 + cone_length))
    vector_end = quaternion.rotate(normal1)

    graphics.cone(molid,
                  tuple(center + vector_start),
                  tuple(center + vector_end),
                  radius=cone_radius, resolution=20)


def get_rotational_matrix(molid, atom_idx1, atom_idx2, axis_start, axis_end, frame_id=-1):
    """
    Get the matrix that rotates the bond between atom 1 and atom 2 onto axis.

    All atoms of selection are rotated on the axis defined by axis_start and
    axis_end. The rotational matrix is then defined by the axis formed of
    atom_idx1 and atom_idx2 (it is rotated so it lies on the former axis).


    Parameters
    ----------
    atom_idx1 : int
        index of first atom

    atom_idx2 : int
        index of second atom

    axis_start : numpy-array
        starting point of the axis (atom_idx1 will be translated onto that point)

    axis_stop : numpy-array
        ending point of the axis; start and end form the axis and therefor
        determine the rotation

    Returns
    -------
    trans_matrix : list {float float ...}
        translational matrix to move atom_idx1 to axis_start
    rot_matrix : list {float float ...}
        rotational matrix to rotate given axis to the desired one

    """
    # get relevant coordinates
    atom_coords1 = vmd_coords_numpy_arrays(molid, frame_id, selection="index {}".format(atom_idx1))
    atom_coords2 = vmd_coords_numpy_arrays(molid, frame_id, selection="index {}".format(atom_idx2))

    # define both axis
    axis1 = atom_coords2 - atom_coords1
    axis1 /= np.linalg.norm(axis1)

    axis2 = axis_end - axis_start

    try:
        axis2 /= np.linalg.norm(axis2)
    except TypeError:
        # vector of length 1 already
        pass

    # angle alpha between axis1 and axis2 (so we can rotate about it)
    alpha = np.arccos(np.dot(axis1, axis2) / (np.linalg.norm(axis1) * np.linalg.norm(axis2)))

    # rotational axis is crossproduct of axis1 and axis2
    rotaxis = np.cross(axis1, axis2)

    # create rotational matrix
    matrix1 = cgt.rotation_matrix(alpha, rotaxis)
    return list(matrix1.flatten())


def measure_geometry(selection, molid, frame=-1, geometry_type="BOND"):
    """
    Measure all bonds in the given atomsel.atomsel object.

    Should be imported from ag_vmd but currently this is not working, therefor
    this duplicate.

    Parameters
    ----------
    selection : atomsel.atomsel
        selection of atoms to measure the bonds for

    Returns
    -------
    bondbond_lengths : dict {tuple{int, int}: float, ...}

    """
    # not sure if this works as intended therefor only imported in this function
    import vmd
    import measure

    geometry_type = geometry_type.upper()

    if geometry_type != "BOND" and geometry_type != "ANGLE" and geometry_type != "DIHEDRAL":
        raise IOError("{} Unknown geometry type!".format(geometry_type))

    geometry_vals = OrderedDict()

    for atom_idx, bond in enumerate(selection.bonds):
        for atom_idx2 in bond:

            # measure bond
            if geometry_type == "BOND":
                # skip reversed bond
                if atom_idx2 == atom_idx:
                    continue

                catms = sorted([atom_idx, atom_idx2])
                cvalue = measure.bond(atom_idx, atom_idx2, molid=molid, frame=frame)[0]
                catms = [i + 1 for i in catms]
                geometry_vals[tuple(catms)] = cvalue
            elif geometry_type == "ANGLE" or geometry_type == "DIHEDRAL":
                for atom_idx3 in selection.bonds[atom_idx2]:
                    # skip reversed bond
                    if atom_idx3 == atom_idx:
                        continue

                    # measure angle
                    if geometry_type == "ANGLE":
                        catms = sorted([atom_idx, atom_idx2, atom_idx3])
                        cvalue = measure.angle(atom_idx, atom_idx2, atom_idx3, molid=molid, frame=frame)[0]
                        catms = [i + 1 for i in catms]
                        geometry_vals[tuple(catms)] = cvalue
                    else:
                        # measure dihedral
                        for atom_idx4 in selection.bonds[atom_idx2]:
                            # skip if just the reversed bond is checked

                            if atom_idx4 == atom_idx3:
                                continue

                            catms = sorted([atom_idx, atom_idx2, atom_idx3, atom_idx4])
                            cvalue = measure.dihedral(atom_idx, atom_idx2, atom_idx3, atom_idx4, molid=molid, frame=frame)[0]
                            catms = [i + 1 for i in catms]
                            geometry_vals[tuple(catms)] = cvalue

    return geometry_vals
