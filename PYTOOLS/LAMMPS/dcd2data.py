#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import ag_unify_md as agum

__version__ = "2017-05-29"

# Convert one or more frames to data-file(s)

parser = argparse.ArgumentParser(prog="dcd2data.py",
                                 formatter_class=argparse.RawTextHelpFormatter,
                                 description="Convert a frame froma a dcd-file" +
                                             "into a data-file.")

parser.add_argument("-d",
                    "--data",
                    metavar="*.lmpdat",
                    required=True,
                    action="store",
                    help="Lammps' data-file(s);\n" + "box-, force field-, " +
                    "topology-parameters must be included!\n" +
                    "http://lammps.sandia.gov/doc/2001/data_format.html",
                    )

parser.add_argument("-dcd",
                    required=True,
                    metavar="*.dcd",
                    action="store",
                    help="Lammps' DCD-file."
                    )

parser.add_argument("-f",
                    "--frame",
                    default=-1,
                    type=int,
                    help="Index (!) of frame to convert (negative indices allowed)."
                    )

parser.add_argument("-o",
                    "--out",
                    default="DEFAULTNAME",
                    action="store",
                    help="Output file name."
                    )

args = parser.parse_args()

frame = args.frame
mydata = agum.Unification()
mydata.read_lmpdat(args.data)
mydata.import_dcd(args.dcd)
mydata.read_frames(frame=None,
                   to_frame=-1,
                   frame_by="index",
                   verbose=False)  # read only the last frame
# convert box to lammps box
print("***Info: Converting box lammps' box format.")
mydata.change_indices(incr=1, mode="increase")
mydata.write_lmpdat(args.out+".lmpdat",
                    frame_id=frame,
                    title="Frame {} of {}".format(frame, args.dcd),
                    cgcmm=True)
