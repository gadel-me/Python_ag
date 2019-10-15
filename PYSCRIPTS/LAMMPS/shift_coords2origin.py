import argparse
import ag_unify_md as agum


if __name__ == '__main__':
    # argument parsing
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("lmpdat",
                        metavar="foo.lmpdat",
                        help="Lammps data file. Provides topology and force field parameters."
                        )

    parser.add_argument("-out",
                        default="foo.lmpdat",
                        metavar="foo.lmpdat",
                        help="Output file."
                        )

    args = parser.parse_args()

    # shift all coordinates so their cog is equal to O (0/0/0)
    mdsys = agum.Unification()
    mdsys.read_lmpdat(args.lmpdat)
    mdsys.transpose_by_cog(copy=False)
    mdsys.change_indices(incr=1, mode="increase")
    mdsys.write_lmpdat(args.out, cgcmm=True)
