from __future__ import print_function, division
import numpy as np
from pyquaternion import Quaternion
import atomsel
import graphics
import color
#import label
import molecule


def draw_angle(molid, frame, atomids, canvas=False, resolution=50, drawcolor="blue"):
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


def label(molid, key, atoms="all", label_color="white", textsize=1.0, offset=(0.0, 0.0, 0.0)):
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


def draw_arrow(molid, start, end, cylinder_radius=0.4, cone_radius=1.0, resolution=50):
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
    graphics.cylinder(molid, tuple(start), tuple(0.9 * end + start),
                      radius=cylinder_radius, resolution=resolution)
    graphics.cone(molid, tuple(0.75 * end + start), tuple(end + start),
                  radius=cone_radius, resolution=resolution)
