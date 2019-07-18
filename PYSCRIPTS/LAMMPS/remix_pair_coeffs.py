#!/usr/bin/env python
from __future__ import print_function, division
import pdb
import argparse
import ag_unify_md as agum

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("lmpt")
    PARSER.add_argument("-epsilon", type=float, default=None)
    PARSER.add_argument("-sigma", type=float, default=None)
    PARSER.add_argument("-atmtype", type=int, default=0)
    PARSER.add_argument("-o", default="foobar.lmpdat")
    ARGS = PARSER.parse_args()

    SYSTEM = agum.Unification()  # force field and topology
    SYSTEM.read_lmpdat(ARGS.lmpt)

    if ARGS.sigma is not None:
        SYSTEM.atm_types[ARGS.atmtype].sigma = ARGS.sigma

    if ARGS.epsilon is not None:
        SYSTEM.atm_types[ARGS.atmtype].epsilon = ARGS.epsilon

    SYSTEM.mix_pair_types(mode="ij")
    SYSTEM.change_indices(incr=1, mode="increase")
    SYSTEM.write_lmpdat(ARGS.o, cgcmm=True)
