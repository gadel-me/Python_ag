#!/usr/bin/env python

import pdb
import argparse
import math
import numpy as np
import md_box as mdb
import ag_unify_md as agum



if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("pwout")
    parser.add_argument("pwin")
    parser.add_argument("-a", nargs=2, type=int, default=(0, 0), metavar=12, help="replicate cell in cell vector a direction this many times.")
    parser.add_argument("-b", nargs=2, type=int, default=(0, 0), metavar=5, help="replicate cell in cell vector b direction this many times.")
    parser.add_argument("-c", nargs=2, type=int, default=(0, 0), metavar=8, help="replicate cell in cell vector c direction this many times.")
    parser.add_argument("-out", default="foobar.pwin", help="Name of output-file.")
    args = parser.parse_args()

    sys  = agum.Unification()
    sys.read_pwin(args.pwin)
    sys.read_pwout(args.pwout)
    sys.replicate_cell(n_start=args.a[0], n_stop=args.a[1], direction="a", frame_id=-1, adjust_box=True)
    sys.replicate_cell(n_start=args.b[0], n_stop=args.b[1], direction="b", frame_id=-1, adjust_box=True)
    sys.replicate_cell(n_start=args.c[0], n_stop=args.c[1], direction="c", frame_id=-1, adjust_box=True)
    sys.write_pwin(-1, args.out)
