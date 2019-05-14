#!/usr/bin/env python
from __future__ import print_function, division
import pdb
import os
import time
import numpy as np
from pyquaternion import Quaternion
from PIL import Image
import ag_unify_md as agum

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


################################################################################
# VMD HELPER FUNCTIONS
################################################################################

def vmd_draw_angle(molid, frame, atomids, canvas=False, resolution=50, drawcolor="blue"):
    #TODO: not working!
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
    id1, id2, id3 = atomids

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
    angle_offset = np.pi/resolution

    for cntr, theta in enumerate(np.arange(0, max_angle, angle_offset)):
        # previous vector
        if cntr == 0:
            vt_pre = vt1 + atm2Crds

        # define quaternion
        q = Quaternion(axis=vtRot, angle=theta)

        # rotate vector
        vt_cur = q.rotate(vt1)
        vt_cur += atm2Crds

        if canvas is True and cntr != 0:
            graphics.triangle(graphics_id, tuple(vt_pre), tuple(atm2Crds), tuple(vt_cur))
        else:
            graphics.cylinder(graphics_id, tuple(vt_pre), tuple(vt_cur), radius=0.01,
                              resolution=20, filled=1)
        vt_pre = vt_cur

    if canvas is True:
        graphics.triangle(graphics_id, tuple(vt_pre), tuple(atm2Crds), tuple(vt2+atm2Crds))
        graphics.material(graphics_id, "Transparent")
    else:
        graphics.cylinder(graphics_id, tuple(vt_pre), tuple(vt2+atm2Crds), radius=0.01,
                          resolution=20, filled=1)


def vmd_label(molid, key, atoms="all", label_color="white", textsize=1.0, offset=(0.0, 0.0, 0.0)):
    """
    Labels atoms by index, charge, etc.

        Input:
           > molid       int; index of molecule to label
           > key         str; attribute to label atom by
           > textsize    float; size of label font
           > offset      (1,3)-list; offset to shift the label by
           > atoms       str or set; "all" or set with atom-ids to label
    """

    # select molecule
    selected_atoms = atomsel.atomsel("all", molid)

    # get its coordinates
    xs = selected_atoms.get("x")
    ys = selected_atoms.get("y")
    zs = selected_atoms.get("z")
    values = selected_atoms.get(key)

    if key == "index":
        values = [i + 1 for i in values]

    # set a label color
    graphics.color(molid, label_color)

    if atoms == "all":
        atoms = range(len(selected_atoms))

    for atom_idx in atoms:
        label_pos = (xs[atom_idx] + offset[0], ys[atom_idx] + offset[1], zs[atom_idx] + offset[2])

        if key == "index":
            labelformat = "{}"
        else:
            labelformat = "{:2.3f}"

        label_text = labelformat.format(values[atom_idx])
        graphics.text(molid, tuple(label_pos), label_text, textsize)


def vmd_draw_arrow(molid, start, end, cylinder_radius=0.4, cone_radius=1.0,
                   cone_length=0.15, resolution=50, double_arrow=False,
                   vmd_material="Basic1Pantone", lcolor="blue"):
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
    """
    graphics.material(molid, vmd_material)
    graphics.color(molid, lcolor)
    p_vector = end - start

    # shift starting point of the vector by cone_length towards the end point
    if double_arrow is True:
        # shorten vector so there is place for the cone at the starting point
        p_start = start + p_vector * cone_length
    else:
        p_start = start

    # shorten vector so there is place for the cone at the ending point
    p_end = start + p_vector * (1 - cone_length)

    graphics.cylinder(molid, tuple(p_start), tuple(p_end), radius=cylinder_radius, resolution=resolution)
    graphics.cone(molid, tuple(p_end), tuple(end), radius=cone_radius, resolution=resolution)
    graphics.cone(molid, tuple(p_start), tuple(start), radius=cone_radius, resolution=resolution)

    # old stuff (working but why?)
    #graphics.cylinder(molid, tuple(start), tuple(0.9 * end + start),
    #                  radius=cylinder_radius, resolution=resolution)
    #graphics.cone(molid, tuple(0.75 * end + start), tuple(end + start), radius=cone_radius, resolution=resolution)


def vmd_render_scene(image_out, image_size=[2000, 2000], renderer="TachyonLOptiXInternal"):
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
    display.set(size=disp_default_size)
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
                      selection="all", scale=0.0, molid=0,
                      vmd_material="Basic1Pantone", style="CPK 1.000000 0.300000 12.000000 12.000000"):
    """
    """
    if molecule.exists(molid) == 0:

        # load molecule in vmd
        if dcd is not None:
            molecule.load(filetype, coordsfile, "dcd", dcd)
        else:
            molecule.load(filetype, coordsfile)

        molrep.modrep(molid, 0, style=style, material=vmd_material, color="Name", sel=selection)
        #trans.scale(scale)
        VMD.evaltcl("scale to {}".format(scale))


def vmd_draw_box(lines, molid=0, vmd_material="Basic1Pantone", lcolor="blue",
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
    graphics.color(molid, lcolor)

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
                          molid=0, vmd_material="Basic1Pantone", textsize=1.0,
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

    origin : tuple
        origin where to place the lattice axis

    molid : int; default = 0
        molid in vmd to draw the lattice axes in

    vmd_material : str
        the material to draw the lattice axes in

    """
    # change vector size
    ucell_a = ucell_a / np.linalg.norm(ucell_a) * 8
    ucell_b = ucell_b / np.linalg.norm(ucell_b) * 8
    ucell_c = ucell_c / np.linalg.norm(ucell_c) * 8

    graphics.material(molid, vmd_material)

    # axis a
    graphics.color(molid, "red")
    vmd_draw_arrow(molid, origin, ucell_a, cylinder_radius=0.4, cone_radius=1.0)
    label_pos = ucell_a * 1.05 + offset + origin
    graphics.text(molid, tuple(label_pos), "a", textsize)

    # axis a
    graphics.color(molid, "green")
    vmd_draw_arrow(molid, origin, ucell_b, cylinder_radius=0.4, cone_radius=1.0)
    label_pos = ucell_b * 1.05 + offset + origin
    graphics.text(molid, tuple(label_pos), "b", textsize)

    # axis c
    graphics.color(molid, "blue")
    vmd_draw_arrow(molid, origin, ucell_c, cylinder_radius=0.4, cone_radius=1.0)
    label_pos = ucell_c * 1.05 + offset + origin
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


def draw_dotted_line(molid, start, end, sradius=0.01, scolor="blue", sdist=0.1, vmd_material="Basic1Pantone"):
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

    scolor : str
        color of the dots

    sdist : float
        distance between points

    """
    graphics.material(molid, vmd_material)
    graphics.color(molid, scolor)

    # vector between start and end
    p_vector = end - start
    length_p_vector = np.linalg.norm(p_vector)

    # normed p_vector
    normed_p_vector = p_vector / length_p_vector

    # number of points
    num_points = int(length_p_vector / sdist)
    print(num_points)

    # first sphere


    for pt in range(num_points):
        csphere = start + normed_p_vector * sdist * pt
        graphics.sphere(molid, center=tuple(csphere), radius=sradius)
