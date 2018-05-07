from __future__ import print_function, division
import math
#import numpy as np
import ag_unify_md as agum
import argparse
import Transformations as cgt
import ag_geometry as agg

# Get angles, dihedrals and impropers for CBZ and write a text file with the given values.
#

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("lmpdat",
                    metavar="foo.lmpdat",
                    help="Lammps data file. Provides topology and force field parameters."
                    )

parser.add_argument("-dcd",
                    help="Lammps dcd-file.")

parser.add_argument("-start",
                    default=0,
                    help="First frame to start analysis with.")

parser.add_argument("-out",
                    default="FOO.txt",
                    help="Name of output file which has angles, dihedrals and impropers.")

args = parser.parse_args()

cbz = agum.Unification()
cbz.read_lmpdat(args.lmpdat)
cbz.import_dcd(args.dcd)
cbz.read_frames(frame=args.start)

with open(args.out, "w") as f_out:
    f_out.write("{:>9}{:>17}{:>17}{:>17}{:>17}\n".format("Step", "ang_C13_N8_C5", "ang_C13_N8_C9", "ang_C9_N8_C5", "ang_b1_N8_b2"))

    for n, frame in enumerate(cbz.ts_coords[args.start:]):

        # omega-1/-2/-3 - angles between carbon and nitrogen
        ang_C13_N8_C5 = cgt.angle_between_vectors(frame[24]-frame[15], frame[24]-frame[7])
        ang_C13_N8_C9 = cgt.angle_between_vectors(frame[24]-frame[15], frame[24]-frame[25])
        ang_C9_N8_C5  = cgt.angle_between_vectors(frame[24]-frame[25], frame[24]-frame[7])

        # phi - angle between benzenes' center of masses
        com_b1 = agg.get_com([frame[14], frame[16], frame[18],
                              frame[20], frame[22], frame[15]],
                             [cbz.atm_types[cbz.atoms[14].atm_key].weigh,
                              cbz.atm_types[cbz.atoms[16].atm_key].weigh,
                              cbz.atm_types[cbz.atoms[18].atm_key].weigh,
                              cbz.atm_types[cbz.atoms[20].atm_key].weigh,
                              cbz.atm_types[cbz.atoms[22].atm_key].weigh,
                              cbz.atm_types[cbz.atoms[15].atm_key].weigh])
        com_b2 = agg.get_com([frame[6], frame[4], frame[0],
                              frame[2], frame[8], frame[7]],
                             [cbz.atm_types[cbz.atoms[6].atm_key].weigh,
                              cbz.atm_types[cbz.atoms[4].atm_key].weigh,
                              cbz.atm_types[cbz.atoms[0].atm_key].weigh,
                              cbz.atm_types[cbz.atoms[2].atm_key].weigh,
                              cbz.atm_types[cbz.atoms[8].atm_key].weigh,
                              cbz.atm_types[cbz.atoms[7].atm_key].weigh])
        ang_b1_N8_b2 = cgt.angle_between_vectors(frame[24]-com_b1, frame[24]-com_b2)

        # improper angle of azepine-n and urea group
        imp_C5_N8_C9_C13 = agg.get_dihedral(frame[7], frame[24], frame[25], frame[15])

        # azepine dihedral
        dih_C13_C12_C15_C7 = agg.get_dihedral(frame[15], frame[14], frame[12], frame[10])
        ang_C13_N8_C5, ang_C13_N8_C9, ang_C9_N8_C5, ang_b1_N8_b2 = [math.degrees(i) for i in (ang_C13_N8_C5, ang_C13_N8_C9, ang_C9_N8_C5, ang_b1_N8_b2)]

        f_out.write("{:>9}{:>17.3f}{:>17.3f}{:>17.3f}{:>17.3f}\n".format(n, ang_C13_N8_C5, ang_C13_N8_C9, ang_C9_N8_C5, ang_b1_N8_b2))
