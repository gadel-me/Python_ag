#!/usr/bin/env python
import pdb
import argparse
import ag_unify_md as agum

parser = argparse.ArgumentParser()

parser.add_argument(
    "pwscf_in",
    metavar="*.pwscf_in",
    help="Quantum Espresso (PW) input file (needed for settings).",
)

parser.add_argument(
    "-pwscf_out",
    default=None,
    metavar="*.pwscf_out",
    help="Quantum Espresso (PW) output file (needed for coordinates and box vectors).",
)

parser.add_argument(
    "-pseudo_dir", default=None, help="directory containing pseudopotential files"
)

parser.add_argument(
    "-restart_mode", default=None, help="Start from scartch or from a previous run."
)

parser.add_argument("-max_seconds", default=None, help="maximum number of seconds")

parser.add_argument(
    "-forc_conv_thr",
    default=None,
    help="Convergence threshold on forces (a.u) for ionic minimization",
)

parser.add_argument(
    "-etot_conv_thr",
    default=None,
    help="Convergence threshold on total energy (a.u) for ionic minimization",
)

parser.add_argument(
    "-ecutwfc",
    type=int,
    default=None,
    metavar="47",
    help="kinetic energy cutoff (Ry) for wavefunctions",
)

parser.add_argument(
    "-ecutrho",
    type=int,
    metavar="323",
    default=None,
    help="Kinetic energy cutoff (Ry) for charge density and potential",
)
parser.add_argument(
    "-conv_thr",
    type=float,
    default=None,
    help="Convergence threshold for selfconsistency: estimated energy error < conv_thr",
)

parser.add_argument(
    "-london_rcut",
    type=int,
    metavar="200",
    default=None,
    help="Kinetic energy cutoff (Ry) for charge density and potential",
)

parser.add_argument(
    "-frame_id",
    type=int,
    default=-1,
    metavar="-1",
    help="ID of Frame from pwscf_out to use coordinates and cell from.",
)

parser.add_argument(
    "-k_points_option",
    type=str,
    default=None,
    metavar="automatic",
    help="K Points option: tpiba | automatic | crystal | gamma | tpiba_b | crystal_b | tpiba_c | crystal_c",
)

parser.add_argument(
    "-k_points_grid",
    type=int,
    nargs=6,
    default=None,
    metavar="-1",
    help="K Points in x, y and z-direction",
)

parser.add_argument(
    "-o",
    default="foo.pwscf_in",
    metavar="*.pwscf_in",
    help="Input file with coordinates and box vectors from 'pwscf_out'",
)

parser.add_argument(
    "-rm_pw_args", default=None, nargs="*", help="Delete pw.x argument",
)


if __name__ == "__main__":
    ARGS = parser.parse_args()
    PW_FILE_HANDLER = agum.Unification()
    PW_FILE_HANDLER.read_pwin(ARGS.pwscf_in)

    # read last frame from pw output file
    if ARGS.pwscf_out is not None:
        PW_FILE_HANDLER.read_pwout(ARGS.pwscf_out)

    if ARGS.pseudo_dir is not None:
        # convert to string with ''
        ARGS.pseudo_dir = f"'{ARGS.pseudo_dir}'"
        PW_FILE_HANDLER.pw_entries["CONTROL"]["pseudo_dir"] = ARGS.pseudo_dir

    if ARGS.max_seconds is not None:
        PW_FILE_HANDLER.pw_entries["CONTROL"]["max_seconds"] = ARGS.max_seconds

    if ARGS.restart_mode == "from_scratch" or ARGS.restart_mode == "restart":
        ARGS.restart_mode = f"'{ARGS.restart_mode}'"
        PW_FILE_HANDLER.pw_entries["CONTROL"]["restart_mode"] = ARGS.restart_mode

    if ARGS.forc_conv_thr is not None:
        PW_FILE_HANDLER.pw_entries["CONTROL"]["forc_conv_thr"] = ARGS.forc_conv_thr

    if ARGS.etot_conv_thr is not None:
        PW_FILE_HANDLER.pw_entries["CONTROL"]["etot_conv_thr"] = ARGS.etot_conv_thr

    # define ecutwfc by user input
    if ARGS.ecutwfc is not None:
        PW_FILE_HANDLER.pw_entries["SYSTEM"]["ecutwfc"] = ARGS.ecutwfc

    # define ecutrho by user input
    if ARGS.ecutrho is not None:
        PW_FILE_HANDLER.pw_entries["SYSTEM"]["ecutrho"] = ARGS.ecutrho

    if ARGS.london_rcut is not None:
        PW_FILE_HANDLER.pw_entries["SYSTEM"]["london_rcut"] = ARGS.london_rcut

    if ARGS.conv_thr is not None:
        PW_FILE_HANDLER.pw_entries["ELECTRONS"]["conv_thr"] = ARGS.conv_thr

    if ARGS.k_points_option is not None:
        PW_FILE_HANDLER.pw_entries["K_POINTS"]["option"] = ARGS.k_points_option

    if ARGS.k_points_grid is not None:
        PW_FILE_HANDLER.pw_entries["K_POINTS"]["k_point_grid"] = ARGS.k_points_grid

    if ARGS.rm_pw_args is not None:
        # list of key pairs to delete
        KEYS_TO_DELETE = []
        # set for faster search
        ARGS.rm_pw_args = set(ARGS.rm_pw_args)

        # get keys and subkeys to delete
        for mainkey, mainvalue in PW_FILE_HANDLER.pw_entries.items():
            for subkey in mainvalue.keys():
                if subkey in ARGS.rm_pw_args:
                    KEYS_TO_DELETE.append((mainkey, subkey))

        for key1, key2 in KEYS_TO_DELETE:
            del PW_FILE_HANDLER.pw_entries[key1][key2]

        del KEYS_TO_DELETE

    PW_FILE_HANDLER.write_pwin(ARGS.frame_id, ARGS.o)
