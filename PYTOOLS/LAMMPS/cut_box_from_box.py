#!/usr/bin/env python
from __future__ import print_function, division
import pdb
import argparse
import numpy as np
import Transformations as cgt
import ag_vectalg as agv
import md_box as mdb
import ag_lammps

#TODO:  Not the expected behavior when shifting the box since the cut_shape
#       function needs an option (like in a dict) to determine if the atoms above (positive distance)
#       or below (negative distance) a plane should be neglected.

"""
Cut a shape by its box vectors from a given lammps file. The shape or the negative
of it may be returned. Cell from lammps data file must be unwrapped.
"""
#===============================================================================
# HELPER FUNCTIONS
#===============================================================================
def cutting(main_sys, cutting_box, inverse=True, frame_id=-1, shift=None):
    """
    Cut the shape of cutting_box from main_sys.

    Parameters
    ----------
    main_sys : ag_unif.Unification()
        system to cut from

    cutting_box : md_box.Box()
        box which defines the cutting shape

    inverse : bool, optional
        inverse the cut

    frame_id : int, optional
        frame index of main_sys to use for cutting
    """
    # convert cutting_box to cartesian
    if cutting_box.boxtype == "lammps":
        cutting_box.box_lmp2cart()

    if shift is None:
        shift = np.array([0, 0, 0])
    elif not isinstance(shift, np.ndarray):
        shift = np.array(shift)
    else:
        pass

    # define box faces
    face_1 = agv.get_plane(cutting_box.crt_a, cutting_box.crt_b, vt_p=shift)  # xy
    face_2 = agv.get_plane(cutting_box.crt_c, cutting_box.crt_a, vt_p=shift)  # xz
    face_3 = agv.get_plane(cutting_box.crt_b, cutting_box.crt_c, vt_p=shift)  # yz
    # multiply with -1 because normal vector must show in the opposite direction
    # the distance between two planes of the same vector is always the third vector (in a parallelepiped at least)
    face_4 = [-1 * i for i in agv.get_plane(cutting_box.crt_a, cutting_box.crt_b, vt_p=cutting_box.crt_c + shift)]
    face_5 = [-1 * i for i in agv.get_plane(cutting_box.crt_c, cutting_box.crt_a, vt_p=cutting_box.crt_b + shift)]
    face_6 = [-1 * i for i in agv.get_plane(cutting_box.crt_b, cutting_box.crt_c, vt_p=cutting_box.crt_a + shift)]
    main_sys.cut_shape(frame_id, inverse, face_1, face_2, face_3, face_4, face_5, face_6)
    cutting_box.box_cart2lmp()


if __name__ == "__main__":
    #===============================================================================
    # ARGUMENT PARSING
    #===============================================================================
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    #-------------------------------------------------------------------------------
    # SYSTEM THAT THE DESIRED SHAPE WILL BE CUT FROM
    #-------------------------------------------------------------------------------
    lmpdat_help = "lammps' data-file which atoms will be cut from."
    dcd_help = "lammps' dcd-file which atoms will be cut from."
    f_help = "Frame index to read."
    rep_help = "replicate cell in cell vector a, b and c direction this many times."
    shft_help = "Translate the coordinates of the (replicated) system to cut from by this translational vector"

    parser.add_argument("lmpdat", metavar="*.lmpdat", help=lmpdat_help)
    parser.add_argument("-dcd", metavar="*.dcd", help=dcd_help)
    parser.add_argument("-f", type=int, default=0, metavar="0", help=f_help)
    parser.add_argument("-rep", default=None, nargs=6, type=int, metavar=4, help=rep_help)
    parser.add_argument("-shft", default=None, nargs=3, type=float, metavar=(1.0, 2.3, 1.2), help=shft_help)

    #-------------------------------------------------------------------------------
    # DEFINE CUTTING SHAPE BY THE BOX OF A GIVEN SYSTEM OR
    # BY THE SYSTEM ITSELF (ATOM COORDINATES AND RADII)
    #-------------------------------------------------------------------------------
    lmpdat_help = "lammps' data-file which defines the cutting shape."
    dcd_help = "lammps' dcd-file which defines the cutting shape."
    f_help = "Frame index of the desired frame"
    rep_help = "replicate cell in cell vector a, b and c direction this many times. Replication takes place before shifting the coordinates."
    shft_help = "Translate the coordinates of the (replicated) system to cut from by this translational vector."
    cut_by_shape_help = "Cut the exact shape of the given system. If a lmpdat_cutshape is given, atomic radii will be considered while cutting."

    parser.add_argument("-lmpdat_cutshape", default=None, metavar="*.lmpdat", help=lmpdat_help)
    parser.add_argument("-dcd_cutshape", default=None, metavar="*.dcd", help=dcd_help)
    parser.add_argument("-f_cutshape", default=-1, type=int, metavar="0", help=f_help)
    parser.add_argument("-rep_cutshape", default=None, nargs=6, type=int, metavar=12, help=rep_help)
    parser.add_argument("-shft_cutshape", default=None, nargs=3, type=float, metavar=("sx, sy, sz"), help=shft_help)
    parser.add_argument("-exact_cutshape", default=False, action="store_true", help=cut_by_shape_help)

    #-------------------------------------------------------------------------------
    # DEFINE CUTTING SHAPE BY THE GIVEN BOX
    #-------------------------------------------------------------------------------
    box_help = "Box vectors in lammps, cartesian or lattice."
    boxtype_help = "Box types: 'lammps' | 'cartesian' | 'lattice'\n" + "Lattice: a, b, c, alpha, beta, gamma in angstrom/degrees\n" + "Lammps: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz"
    parser.add_argument("-box", default=None, type=float, nargs="*", help=box_help)
    parser.add_argument("-boxtype", help=boxtype_help)

    #-------------------------------------------------------------------------------
    # DEFINE CUTTING SHAPE BY PLANES
    #-------------------------------------------------------------------------------
    cut_planes_help = "Planes that define the box which will be cut - TBD."
    parser.add_argument("-cut_planes", type=float, nargs="*", help=cut_planes_help)

    #-------------------------------------------------------------------------------
    # CUTTING OPTIONS
    #-------------------------------------------------------------------------------
    inverse_help = "Selection: box (normal), negative of box (inverted)"
    shft_cut_result = "Shift the cut out coordinates by these 3 vectors (index 0, 1, 2: a; 3, 4, 5: b; 6, 7, 8: c"
    enlarge_cut_help = "Enlarge the box of the cuttings shape by this vector"
    parser.add_argument("-inverse", default=False, action="store_true", help=inverse_help)
    parser.add_argument("-shft_cut_result", nargs=3, type=float, default=(0.0, 0.0, 0.0), metavar=("sx, sy, sz"), help=shft_cut_result)
    parser.add_argument("-scale_cut", default=None, nargs=9, type=float, metavar=2.42, help=enlarge_cut_help)

    parser.add_argument("-out", default="test.lmpdat")

    args = parser.parse_args()

    #===============================================================================
    # PREPARE THE MAIN SYSTEM
    #===============================================================================
    sys_cutfrom = ag_lammps.read_lmpdat(args.lmpdat, dcd=args.dcd, frame_idx_start=args.f, frame_idx_stop=args.f)

    # replicate cell if desired
    if args.rep_cutshape is not None:
        sys_cutfrom.replicate_cell(n_start=args.rep[0], n_stop=args.rep[1], direction="a", frame_id=-1, adjust_box=True)
        sys_cutfrom.replicate_cell(n_start=args.rep[2], n_stop=args.rep[3], direction="b", frame_id=-1, adjust_box=True)
        sys_cutfrom.replicate_cell(n_start=args.rep[4], n_stop=args.rep[5], direction="c", frame_id=-1, adjust_box=True)
        sys_cutfrom.fetch_molecules_by_bonds()
        sys_cutfrom.mols_to_grps()

    # shift sys_cutfrom by given vector
    if args.shft is not None:
        args.shft = np.array(args.shft)
        M_shft = cgt.translation_matrix(args.shft)
        atm_idxs = range(len(sys_cutfrom.atoms))
        sys_cutfrom.mm_atm_coords(-1, M_shft, False, *atm_idxs)
        del (atm_idxs, M_shft)

    #===============================================================================
    # PREPARE THE CUTTING SHAPE SYSTEM
    #===============================================================================
    if args.lmpdat_cutshape is not None or args.dcd_cutshape is not None:
        sys_cutshape = ag_lammps.read_lmpdat(args.lmpdat_cutshape, dcd=args.dcd_cutshape, frame_idx_start=args.f_cutshape, frame_idx_stop=args.f_cutshape)

        # replicate cell if desired
        if args.rep_cutshape is not None:
            sys_cutshape.replicate_cell(n_start=args.rep_cutshape[0], n_stop=args.rep_cutshape[1], direction="a", frame_id=-1, adjust_box=True)
            sys_cutshape.replicate_cell(n_start=args.rep_cutshape[2], n_stop=args.rep_cutshape[3], direction="b", frame_id=-1, adjust_box=True)
            sys_cutshape.replicate_cell(n_start=args.rep_cutshape[4], n_stop=args.rep_cutshape[5], direction="c", frame_id=-1, adjust_box=True)

        # shift the coordinates (not the box!) of sys_cutshape by given vector
        # (for cutting the exact shape of sys_cutshape)
        if args.shft_cutshape is not None:
            args.shft_cutshape = np.array(args.shft_cutshape)
            M_shft_cut = cgt.translation_matrix(args.shft_cutshape)
            atm_idxs = range(len(sys_cutshape.ts_coords[-1]))
            sys_cutshape.mm_atm_coords(-1, M_shft_cut, False, *atm_idxs)
            del (atm_idxs, M_shft_cut)

    else:
        sys_cutshape = None

    #===============================================================================
    # PREPARE BOX WHICH SOLELY DEFINES THE SHAPE
    #===============================================================================
    if args.box is not None:
        import math

        if args.boxtype == "lammps":
            box = mdb.Box(boxtype=args.boxtype, lmp_xlo=args.box[0], lmp_xhi=args.box[1], lmp_ylo=args.box[2], lmp_yhi=args.box[3], lmp_zlo=args.box[4], lmp_zhi=args.box[5], lmp_xy=args.box[6], lmp_xz=args.box[7], lmp_yz=args.box[8])
            #box.box_lmp2cart()  # change box vectors to lattice
        elif args.boxtype == "lattice":
            box = mdb.Box(boxtype=args.boxtype, ltc_a=args.box[0], ltc_b=args.box[1], ltc_c=args.box[2], ltc_alpha=math.radians(args.box[3]), ltc_beta=math.radians(args.box[4]), ltc_gamma=math.radians(args.box[5]))
            box.box_lat2cart()
        else:
            box = mdb.Box(boxtype=args.boxtype, crt_a=[args.box[0], args.box[1], args.box[2]], crt_b=[args.box[3], args.box[4], args.box[5]], crt_c=[args.box[6], args.box[7], args.box[8]])
    else:
        box = None

    # equalize box variable
    if sys_cutshape is not None:
        box = sys_cutshape.ts_boxes[-1]

    if box.boxtype == "lammps":
        box.box_lmp2cart()

    #===============================================================================
    # CUT THE BOX
    #===============================================================================
    # change box by vector
    if args.scale_cut is not None:
        box.crt_a += np.array(args.scale_cut[0:3])
        box.crt_b += np.array(args.scale_cut[3:6])
        box.crt_c += np.array(args.scale_cut[6:9])

    cutting(sys_cutfrom, box, args.inverse, -1, shift=args.shft_cutshape)

    #===============================================================================
    # CUT THE BOX
    #===============================================================================
    # assign new box vectors to the existing box
    if args.inverse is True:
        sys_cutfrom.ts_boxes[-1] = box

    #pdb.set_trace()
    sys_cutfrom.refresh()
    sys_cutfrom.change_indices(incr=1, mode="increase")

    if args.shft_cut_result is not None:
        M_shft_cut = cgt.translation_matrix(args.shft_cut_result)
        atm_idxs = range(len(sys_cutfrom.ts_coords[-1]))
        sys_cutfrom.mm_atm_coords(-1, M_shft_cut, False, *atm_idxs)
        del (atm_idxs, M_shft_cut)

    print("Writing lammps data file...")
    sys_cutfrom.write_lmpdat(args.out, frame_id=-1, title="Cut box", cgcmm=True)
