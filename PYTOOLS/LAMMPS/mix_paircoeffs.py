import ag_unify_md as agum
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("dat")

args = parser.parse_args()

# Read Data --------------------------------------------------------------------
mydata = agum.Unification()
mydata.read_lmpdat(args.dat)
mydata.mix_pair_types(mode="ij", to_file="{}.lmpcfg".format(args.dat))
