#!/usr/bin/env python
from __future__ import print_function, division
import readline  # necessary for raw_input and using arrow keys
import argparse
import numpy as np
import ag_unify_log as agul
import plot_clmsv_helper as pcsvh


parser = argparse.ArgumentParser(prog="ag_plot_lmplog.py",
                                 formatter_class=argparse.RawTextHelpFormatter,
                                 description="Plot data from log.lammps.")
group_dor2 = parser.add_mutually_exclusive_group(required=False)

parser.add_argument("clmsv",
                    metavar="foo.clmsv",
                    action="store",
                    help="Column separated files containing values " +
                         "for several keywords"
                    )


parser.add_argument("-xkey",
                    "--thermo_x",
                    dest="xkey",
                    metavar="steps|press|temp",
                    action="store",
                    default="Step",
                    required=False,
                    help="output(s) from thermo-setup as x-value."
                    )

parser.add_argument("-ykey",
                    "--thermo_y",
                    dest="ykey",
                    metavar="Step|Press|Temp, etc.",
                    action="store",
                    default=None,
                    required=False,
                    help="output(s) from thermo-setup as y-value."
                    )

parser.add_argument("-fst",
                    action="store",
                    type=int,
                    default=0,
                    help="Thermo output from first frame."
                    )

parser.add_argument("-lst",
                    action="store",
                    type=int,
                    default=-1,
                    help="Thermo output from last frame."
                    )

group_dor2.add_argument("-linreg",
                        "--linear_regression",
                        dest="linreg",
                        action="store_true",
                        default=False,
                        help="plot the linear regression"
                        )

group_dor2.add_argument("-histo",
                        "--histogram",
                        dest="histo",
                        action="store_true",
                        default=False,
                        help="plot data as histogram"
                        )

parser.add_argument("-fit",
                    "--fit_function",
                    dest="fit",
                    help="use a function to fit plotted data."
                    )

parser.add_argument("-mean",
                    action="store_true",
                    default=False,
                    help="calculate the mean from '-start' to '-stop' of the y-values"
                    )

parser.add_argument("-median",
                    action="store_true",
                    default=False,
                    help="calculate the median from '-start' to '-stop' of the y-values"
                    )

args = parser.parse_args()

# read clmsv
clmsv = agul.LogUnification()
clmsv.read_clmsv(args.clmsv)
print("***Warning: Currently only the first entry of each clmsv will be plotted!")
data = clmsv.data[0]
keys = data.keys()

# /// plot data ///
keep_plotting = "y"
first_time = True

previous_xkey = args.xkey
previous_ykey = args.ykey
previous_fst  = args.fst
previous_lst  = args.lst
num_frames = len(data[args.xkey])

while keep_plotting in ["Yes", "Y", "yes", "y", ""]:
    # switches to check for right input
    wrong_fst  = True
    wrong_lst  = True
    wrong_xkey = True
    wrong_ykey = True

    print("Thermodynamic data available for {} frames.".format(num_frames))

    # 1st frame
    while wrong_fst is True:
        args.fst = pcsvh.ask4frame("first", num_frames)

        if args.fst == "":
            args.fst = previous_fst

        args.fst  = int(args.fst)

        if args.fst < -1 or args.fst > num_frames:
            print("{} is not a valid frame number!".format(args.fst))
        else:
            print("First frame: {}.\n".format(args.fst))
            break

    # last frame
    while wrong_lst is True:
        args.lst = pcsvh.ask4frame("last", num_frames)

        if args.lst == "":
            args.lst = previous_lst

        args.lst = int(args.lst)

        if args.lst > num_frames:
            print("{} beyond last frame!".format(args.lst))
        elif args.lst < args.fst and args.lst != -1:
            print("Last frame ({}) smaller than first frame! ({})".format(args.lst, args.fst))
        else:
            print("Last frame: {}.\n".format(args.lst))
            break

    # x-key
    while wrong_xkey is True:
        args.xkey = pcsvh.ask4keyword("x-values", keys)
        if args.xkey == "":
            args.xkey = previous_xkey

        if args.xkey not in keys:
            print("{} not a valid keyword!\n".format(args.xkey))
        else:
            print("\nUsing values from {} (x).".format(args.xkey))
            break

    xvals = data[args.xkey][args.fst:args.lst]

    # y-key
    if args.histo is False:
        while wrong_ykey:
            args.ykey = pcsvh.ask4keyword("y-values", keys)

            if args.ykey == "":
                args.ykey = previous_ykey

            if args.ykey not in keys:
                print("{} not a valid keyword!\n".format(args.ykey))
            else:
                print("\nUsing values from {} (y).".format(args.ykey))
                break

        yvals = data[args.ykey][args.fst:args.lst]

    # plotting
    if args.histo is True:
        xmean = np.mean(xvals)
        print("Mean of {}: {}".format(args.xkey, xmean))

        print("Plotting: {} for steps {} to {}".format(args.xkey, args.fst,
                                                       args.lst))
        pcsvh.plot_histogram(xvals, args.xkey)
    else:
        pcsvh.plot_xy(xvals, yvals, args.xkey, args.ykey)

    previous_xkey = args.xkey
    previous_ykey = args.ykey
    previous_fst  = args.fst
    previous_lst  = args.lst
    # continue plotting?
    keep_plotting = raw_input("Plot something else? Y(es)/n(o): > ")

    if keep_plotting in ["Yes", "Y", "yes", "y", ""]:
        print("\nNext plot...")
    else:
        print("Closing...")
        break
