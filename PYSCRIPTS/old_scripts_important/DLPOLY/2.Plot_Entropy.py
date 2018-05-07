#!/home/angad/Programs/anaconda/bin/python

from __future__ import print_function, division
import sys
sys.path.append("/home/gadelmeier/Python/myPYMODULES/DLPOPLY_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/OTHER_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/MATH_MODULES/")
import numpy as np
#import scipy.stats as stats
import DL_HISTORY_class8_new_approaches_8_2 as dlhc
#import my_lin_alg_compendium as mlc
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
import tqdm
import collections

# USAGE:
# this_script, RESULTFILE1, FIELDFILE1, RESULTFILE2, FIELDFILE2, ...


def plot_histogram(data_dictionary, label=None, plot_to_png=False, debug=False):
    # PREPARE BOUNDARY-CONDITIONS FOR PLOTTING
    # --> x-range, y-range, number of subplots, number of bins
    all_x_min = 9e20  # dummy value for minimum of all data lists
    all_x_max = -9e20  # dummy value for maximum of all data lists
    most_data_pts = None
    total_num_data_plots = 0  # number of plots

    for cur_key in data_dictionary:
        data = data_dictionary[cur_key]

        # skip empty data arrays
        if data:
            pass
        else:
            continue

        # increase number plots
        total_num_data_plots += 1
        # find lowest and biggest value
        cur_min_x = min(data)
        cur_max_x = max(data)

        if debug:
            print(cur_min_x, cur_max_x)

        if cur_min_x < all_x_min:
            all_x_min = cur_min_x

            if debug:
                print("New minium found", all_x_min)

        if cur_max_x > all_x_max:
            all_x_max = cur_max_x

            if debug:
                print("New maximum found", all_x_max)

        cur_data_points = len(data)

        if cur_data_points > most_data_pts:
            most_data_pts = cur_data_points

    # SET NUMBER OF BINS ///////////////////////////////////////////////////////
    sqrt_bins = math.ceil(
        np.sqrt(most_data_pts)
    )

    # reduce number of bins if too many
    if sqrt_bins > 200:
        sqrt_bins = 200

    if debug:
        print("Number of most data points: ", most_data_pts)
        print("Square rooted bins: ", sqrt_bins)

    # get evenly wide bins +++++++++++++++++++++++++++++++++++++++++++++++++++++
    x_range = np.linspace(all_x_min, all_x_max, sqrt_bins)

    # get min and max y-values +++++++++++++++++++++++++++++++++++++++++++++++++
    all_y_min = 9e20
    all_y_max = -9e20  # dummy value for maximum of all y-values

    for dummy_key2 in data_dictionary:
        dummy_data2 = data_dictionary[dummy_key2]

        if dummy_data2:
            pass
        else:
            continue  # exclude empty bins

        dummy_n, dummy_bins = np.histogram(dummy_data2, bins=x_range,
                                           normed=True)

        cur_min_hist_y = min(dummy_n+dummy_n/10)
        cur_max_hist_y = max(dummy_n+dummy_n/10)

        if cur_min_hist_y < all_y_min:
            all_y_min = cur_min_hist_y
        if cur_max_hist_y > all_y_max:
            all_y_max = cur_max_hist_y

    # PROCESS RAW DATA AND PLOT ////////////////////////////////////////////////
    fig = plt.figure()
    cur_plot_nr = 0  # order rank of subplots

    # Iterate through dictionary
    for key in data_dictionary:
        data = data_dictionary[key]

        if data:
            pass
        else:
            continue  # data empty? process next data container

        # add subplots
        cur_plot_nr += 1
        # rows, columns, order rank
        fig.add_subplot(total_num_data_plots, 1, cur_plot_nr)

        # DEFINE CURRENT SUBPLOT ///////////////////////////////////////////////

        # define frame conditions for plotting
        plt.xlim(all_x_min, all_x_max)
        plt.ylim(all_y_min, all_y_max)
        plt.xlabel("delta", fontsize=10, fontweight='bold')
        plt.ylabel("Cumulative probability/\n or probability density",
                   fontsize=10, fontweight='bold')
        plt.grid(True)  # Set a grid

        # get some statistics data for pdf and histogram
        cur_mu = np.mean(data)
        cur_sigma = np.std(data)
        #cur_median = np.median(data)

        # plot histogram
        n, bins, patches = plt.hist(data, bins=x_range,
                                    normed=True, cumulative=False,
                                    histtype="step")

        if debug:
            print("Max n: ", max(n))
            if bins.all() == x_range.all():
                print("bins and x_range are equal!")
            else:
                print("bins and x_range are not equal!")

        # plot best fit line for given x_range, mu and sigma
        pdf_x_values = bins
        pdf_y_values = mlab.normpdf(pdf_x_values, cur_mu, cur_sigma)
        plt.plot(pdf_x_values, pdf_y_values, "r--", linewidth=0.5)

        # line that goes through mu (=max. of gauss-function)
        plt.axvline(x=cur_mu, linewidth=0.5, color="red", alpha=1.0)
        # line that goes through median
        #plt.axvline(x=cur_median, linewidth=0.5, color="orange", alpha=1.0)

        # Set a title for current subplot
        plt.title(
            "Histogram of %r - %s  $\mu=%.4f$ $\sigma=%.4f$" % (key, label, cur_mu, cur_sigma),
            fontweight='bold', fontsize=11
        )
        # ?
        plt.subplots_adjust(left=0.15)

    # Show Plot
    plt.show()
    # Save Plot to PNG-File
    if plot_to_png:
        fig.savefig("%s.png" % (label), dpi=200, bbox_inches="tight")


# FETCH DATA FROM FILES ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
result_files = sys.argv[1::2]  # every 2nd value, starting with second item in list
field_files = sys.argv[2::2]  # every 2nd value, starting with third item in list
# File progression
num_result_files = range(len(result_files))
pbar = tqdm.tqdm(num_result_files)

moltypes_com = {}
moltypes_chi = {}
moltypes_phi = {}
moltypes_psi = {}

for resultfile, fieldfile, tqdmstep in zip(result_files, field_files, pbar):
    # READ FIELD FILES /////////////////////////////////////////////////////////
    finfo = dlhc.FieldFile(fieldfile)
    finfo.read_field()

    # PREPARE CONTAINERS ///////////////////////////////////////////////////////
    moltypes = [i.molname for i in finfo.molecule_types]  # get molnames
    for moltype in finfo.molecule_types:
        # molecule type names must vary in each FIELD File
        # e.g. if "TIP3P Water" is a common molname in all FIELD Files, but
        # data for "TIP3P Water" will not be separated but appended and cannot
        # be distinguished in the end
        moltypes_com[moltype.molname] = []
        moltypes_chi[moltype.molname] = []
        moltypes_phi[moltype.molname] = []
        moltypes_psi[moltype.molname] = []

    # READ RESULT FILES ////////////////////////////////////////////////////////
    cur_moltype = None
    with open(resultfile, "r") as result_file:
        for line in result_file:
            line = line.rstrip()

            if line in moltypes:
                cur_moltype = line

            elif "COM" in line:  # omit plain text line in RESULT File
                pass
            else:
                line = line.split()
                moltypes_com[cur_moltype].append(
                    float(line[0])
                )
                moltypes_chi[cur_moltype].append(
                    float(line[1])
                )
                moltypes_phi[cur_moltype].append(
                    float(line[2])
                )
                moltypes_psi[cur_moltype].append(
                    float(line[3])
                )

    pbar.set_description("Processing %s" % resultfile)
    pbar.update(1)

# Order dicts by keys ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
moltypes_com = collections.OrderedDict(
    sorted(
        moltypes_com.items()
    )
)
moltypes_chi = collections.OrderedDict(
    sorted(
        moltypes_chi.items()
    )
)
moltypes_phi = collections.OrderedDict(
    sorted(
        moltypes_phi.items()
    )
)
moltypes_psi = collections.OrderedDict(
    sorted(
        moltypes_psi.items()
    )
)

# Plot /////////////////////////////////////////////////////////////////////////
plot_histogram(moltypes_com, label="COM", plot_to_png=True)
plot_histogram(moltypes_chi, label="CHI", plot_to_png=True)
plot_histogram(moltypes_phi, label="PHI", plot_to_png=True)
plot_histogram(moltypes_psi, label="PSI", plot_to_png=True)
