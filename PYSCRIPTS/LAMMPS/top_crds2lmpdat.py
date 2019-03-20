#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import math
import numpy as np
import md_box as mdb
import ag_unify_md as agum
import pdb

"""
Input is as follows:

    Coordinates: > Amber (inpcrd, mdcrd)
                 > Lammps (data-, molecule-file
                 > xyz-file\n\
    Topology and ForceField: > Amber (prmtop)
                             > Lammps (data-file)

Caveat: Order of atoms in topology and coordinates must be the same!
If n-times of topology atoms are in the coordinates file, the topology will be
duplicated by factor n (n is int).
"""


def get_cell_from_cif(cif_file):
    """Get lattice entities from a cif file."""
    def get_number_from_line(input_line):
        """Get lattice entities from line."""
        input_line = input_line.split()
        number = float(input_line[1].split("(")[0])
        return number

    with open(cif_file) as f_open:
        line = f_open.readline()

        while line != "":
            if "_cell_length_a" in line:
                cif_a = get_number_from_line(line)
            elif "_cell_length_b" in line:
                cif_b = get_number_from_line(line)
            elif "_cell_length_c" in line:
                cif_c = get_number_from_line(line)
            elif "_cell_angle_alpha" in line:
                cif_alpha = get_number_from_line(line)
                cif_alpha = np.radians(cif_alpha)
            elif "_cell_angle_beta" in line:
                cif_beta = get_number_from_line(line)
                cif_beta = np.radians(cif_beta)
            elif "_cell_angle_gamma" in line:
                cif_gamma = get_number_from_line(line)
                cif_gamma = np.radians(cif_gamma)
            else:
                pass

            line = f_open.readline()

    box = mdb.Box(boxtype="lattice", ltc_a=cif_a, ltc_b=cif_b, ltc_c=cif_c,
                  ltc_alpha=cif_alpha, ltc_beta=cif_beta, ltc_gamma=cif_gamma)
    return box


def get_other_stuff_from_cif(cif_file):
    """Get lattice entities from a cif file."""
    def get_number_from_line(input_line):
        """Get lattice entities from line."""
        input_line = input_line.split()
        number = float(input_line[1].split("(")[0])
        return number

    with open(cif_file) as f_open:
        line = f_open.readline()

        while line != "":
            if "_cell_volume" in line:
                cell_volume = get_number_from_line(line)
            elif "_cell_formula_units_Z" in line:
                cell_formula_units = get_number_from_line(line)
            elif "_cell_measurement_temperature" in line or "_diffrn_ambient_temperature" in line:
                cell_measurement_temperature = get_number_from_line(line)
            elif "_exptl_crystal_density_diffrn" in line:
                density = get_number_from_line(line)
            else:
                pass

            line = f_open.readline()

        return {"cell_volume": cell_volume, "cell_formula_units": cell_formula_units, "cell_measurement_temperature": cell_measurement_temperature, "density": density}


if __name__ == "__main__":
    # Argument Parsing -------------------------------------------------------------
    parser = argparse.ArgumentParser(
        prog="top_crds2lmpdat.py",
        formatter_class=argparse.RawTextHelpFormatter,
        description="Merge topology and coordinates from different files." +
                    "Replicates topology information automatically if necessary." +
                    "Same order in coordinates file and topology/force field file is crucial."
    )

    # args for topology/forcefield and coordinates
    topo_ff = parser.add_mutually_exclusive_group(required=True)
    coords = parser.add_mutually_exclusive_group(required=True)
    cell = parser.add_mutually_exclusive_group()

    topo_ff.add_argument("-prmtop",
                         default=None,
                         metavar="foo.prmtop",
                         help="Amber file. Provides topology and force field parameters.",
                         )

    topo_ff.add_argument("-lmpt",
                         "--lmpdat_top",
                         default=None,
                         metavar="foo.lmpdat",
                         help="Lammps data file. Provides topology and force field parameters."
                         )

    coords.add_argument("-xyz",
                        default=None,
                        metavar="bar.xyz",
                        help="XYZ file. Provides coordinates.",
                        )

    coords.add_argument("-pwout",
                        default=None,
                        metavar="bar.pw_out",
                        help="PW-output file. Provides coordinates.",
                        )

    coords.add_argument("-inpcrd",
                        default=None,
                        metavar="bar.inpcrd",
                        help="Amber file. Provides coordinates.")

    coords.add_argument("-lmpc",
                        "--lmpdat_coords",
                        default=None,
                        metavar="bar.lmpdat",
                        help="Lammps data file. Provides coordinates" +
                             "(force field and topology information will be ignored)."
                        )

    coords.add_argument("-gau_in",
                        default=None,
                        metavar="foo.gau",
                        help="Gaussian input file with coordinates.")

    parser.add_argument("-sysname",
                        metavar="UNK",
                        default="UNK",
                        help="Name of the molecule/system")

    cell.add_argument("-cell",
                      nargs=6,
                      metavar=("a", "b", "c", "alpha", "beta", "gamma"),
                      type=float,
                      default=None,
                      help="Box vectors a, b, c and box angles alpha, beta, gamma" +
                           "(in degrees) in that particular order."
                      )

    cell.add_argument("-cbg",
                      "--cell_by_guess",
                      action="store_true",
                      help="Guess cell by atomic coordinates."
                      )

    cell.add_argument("-cbc",
                      "--cell_by_coordfile",
                      action="store_true",
                      help="Get cell information from coordinate file (currently only reads cell from lammps data file)."
                      )

    cell.add_argument("-cbcif",
                      "--cell_by_ciffile",
                      default=None,
                      help="Get cell information from a cif file")

    parser.add_argument("-out",
                        metavar="foobar.lmpdat",
                        default="foobar.lmpdat",
                        help="Name of output-file."
                        )

    args = parser.parse_args()

    # Modeling ---------------------------------------------------------------------
    sys_coords = agum.Unification()  # coordinates
    sys_topo_ff = agum.Unification()  # force field and topology

    # read topology and force field parameters
    if args.prmtop is not None:
        sys_topo_ff.read_prmtop(args.prmtop)
        sys_topo_ff.ui_convert_units(energy_unit_out='eV',
                                     ang_unit_out="deg",
                                     cvff_style=False)
        sys_topo_ff.mix_pair_types(mode="ij")
    elif args.lmpdat_top is not None:
        sys_topo_ff.read_lmpdat(args.lmpdat_top, energy_unit="eV", angle_unit="deg")
    else:
        raise Warning("This should not have happened (1).")

    # read coordinates
    if args.lmpdat_coords is not None:
        sys_coords.read_lmpdat(args.lmpdat_coords)
    elif args.xyz is not None:
        sys_coords.read_xyz(args.xyz)
    elif args.inpcrd is not None:
        sys_coords.read_inpcrd(args.inpcrd)
    elif args.pwout is not None:
        sys_coords.read_pwout(args.pwout)
    elif args.gua_in is not None:
        sys_coords.read_gau(args.gau_in)
    else:
        raise Warning("This should not have happened (2).")

    # replace or append coordinates if they have already been provided/ not provided
    try:
        sys_topo_ff.ts_coords[0] = sys_coords.ts_coords[0]
    except IndexError:
        sys_topo_ff.ts_coords.append(sys_coords.ts_coords[0])

    # replicate topology if necessary
    num_topo_atms   = len(sys_topo_ff.atoms)
    num_topo_coords = len(sys_topo_ff.ts_coords[0])
    n = num_topo_coords/num_topo_atms

    if n.is_integer():
        # replicate topology
        sys_topo_ff.add_topology_replicate(int(n-1), refresh_bonds=True)
    else:
        raise Warning("Number of atoms in coordinates file is not a multiple of atoms in topology file!")

    # box --------------------------------------------------------------------------
    # check if coordinates file has a cell given
    if args.cell_by_coordfile is True:
        if sys_coords.ts_boxes == []:
            print("***Warning: Coordinate file has no box info! Guessing cell by coordinates")
            args.cell_by_guess = True
            args.cell_by_coordfile = False

    # cell coordinates by user
    if args.cell is not None:
        a, b, c = args.cell[0:3]
        alpha, beta, gamma = [math.radians(i) for i in args.cell[3:]]

        # replace box if given or append if none was given before
        if sys_topo_ff.ts_boxes == []:
            sys_topo_ff.ts_boxes.append(mdb.Box(boxtype="lattice",
                                                ltc_a=a, ltc_b=b, ltc_c=c,
                                                ltc_alpha=alpha, ltc_beta=beta,
                                                ltc_gamma=gamma))
        else:
            sys_topo_ff.ts_boxes[0] = mdb.Box(boxtype="lattice",
                                              ltc_a=a, ltc_b=b, ltc_c=c,
                                              ltc_alpha=alpha, ltc_beta=beta,
                                              ltc_gamma=gamma)

        sys_topo_ff.ts_boxes[0].box_lat2lmp()
    # guess cell size if none was defined or user wants it that way explicitly
    elif args.cell_by_guess is True or sys_topo_ff.ts_boxes == []:
        sys_topo_ff.def_boxes_by_coords()
    elif args.cell_by_coordfile is True:
        sys_topo_ff.ts_boxes[0] = sys_coords.ts_boxes[0]
    elif args.cell_by_ciffile is not None:
        sys_topo_ff.ts_boxes[0] = get_cell_from_cif(args.cell_by_ciffile)
    else:
        pass  # keep given cell

    # group atoms by bonds
    sys_topo_ff.fetch_molecules_by_bonds()
    sys_topo_ff.mols_to_grps()
    # change indices to start with 1 (everything is newly ordered,
    # internally starting with 0 from 1st atom given in coordinates)
    sys_topo_ff.change_indices(incr=1, mode="increase")
    sys_topo_ff.write_lmpdat(args.out, frame_id=0, title=args.sysname, cgcmm=True)

    # just for debugging
    # already done by write lammps data routine - just for testing
    #sys_topo_ff.ts_boxes[0].box_lat2lmp()
    #sys_topo_ff.ts_boxes[0].box_lmp2lat()
    #print(sys_topo_ff.ts_boxes[0].__dict__)
    #print(np.degrees(sys_topo_ff.ts_boxes[0].ltc_alpha))
    #print(np.degrees(sys_topo_ff.ts_boxes[0].ltc_beta))
    #print(np.degrees(sys_topo_ff.ts_boxes[0].ltc_gamma))
