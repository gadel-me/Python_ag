import re
import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse
import ag_statistics as ags
from restrain_scan import norm_energy
from ag_statistics import gnuplot_gaussfit
from pantone_colors import pantone_colors_hex
import pdb


def read_normed_output(filename):
    """
    Read a two column file with one header line.

    All lines containing # will be skipped.

    Parameters
    ----------
    filename : str

    Returns
    -------
    entities : list
        first column
    energies : list
        second column

    """
    entities = []
    energies = []

    with open(filename) as open_file:
        line = open_file.readline()

        while line != '':
            line = open_file.readline()

            try:
                entity, energy = [float(number) for number in line.split()[:2]]
                entities.append(entity)
                energies.append(energy)
            except ValueError:
                pass

    return (entities, energies)


def plot_all(ref_file, data_files, xlabel=None, title=None, substract=False, x_offset=None):
    """
    Plot all data.

    Parameters
    ----------
    ref_file : str
    data_files : list of str
    xlabel : None or str

    """
    fig = plt.figure()

    def plot_data_files():
        """
        Plot the given data and calculating chi square.

        Calculate chi square error and extrapolate the data to the reference.
        """
        # colormap
        norm = matplotlib.colors.Normalize(vmin=0, vmax=len(data_files))
        cmap = matplotlib.cm.get_cmap('Dark2')
        cmap_reds = matplotlib.cm.get_cmap('autumn')

        for index, data_file in enumerate(data_files):
            iteration = int(re.findall(r"\d+", data_file)[-1])

            if "without_entity" in data_file:
                addname = "w/o entity"
            else:
                addname = ""

            data_x, data_y = read_normed_output(data_file)

            # only for dihedrals
            #data_x = [i + 180 if i < 0 else i - 180 for i in data_x]
            #data_x = [i -4 for i in data_x]

            interp_y = np.interp(ref_x, data_x, data_y)
            chi_square_error = ags.chi_square_error(ref_y, interp_y)
            label = r"{:> 6} {} - $\chi^2$-err:{:> 8.7f} eV".format(iteration, addname, chi_square_error)
            color = cmap(norm(index))
            if x_offset is not None:
                plt.plot([i + x_offset if i < 0 else i for i in data_x], data_y, linestyle=lstyle, marker=marker, linewidth=0.2, markersize=0.5, color=color)
                plt.plot([i + x_offset if i < 0 else i for i in ref_x], interp_y, linestyle=lstyle, marker=marker, linewidth=0.2, color=color, label=label)
            else:
                plt.plot(data_x, data_y, linestyle=lstyle, marker=marker, linewidth=0.2, markersize=0.5, color=color)
                plt.plot(ref_x, interp_y, linestyle=lstyle, marker=marker, linewidth=1.0, color=color, label=label)

            # substract the reference from the current iteration
            if substract is True:
                red = cmap_reds(norm(index))
                data_delta = [abs(i) - abs(j) for i, j in zip(ref_y, interp_y)]

                if x_offset is not None:
                    plt.plot([i + x_offset if i < 0 else i for i in ref_x], data_delta, linestyle="--", linewidth=0.2, marker=marker, color=red, label=r"{:> 6} vs. ab initio: $\Delta$ {}".format(iteration, addname))
                else:
                    #plt.plot(ref_x, data_delta, linestyle="--", linewidth=0.2, marker=marker, color=red, label=r"{:> 6} vs. ab initio: $\Delta$ {}".format(iteration, addname))
                    plt.plot(ref_x, data_delta, linestyle="--", linewidth=1.0, marker=marker, color=red, label=r"{:> 6} vs. ab initio: $\Delta$ {}".format(iteration, addname))

                # write delta energies to a file for further plotting
                print("***Info: Writing delta file!")
                with open(data_file.rstrip(".txt") + "_delta.txt", "w") as fin:
                    fin.write("     Angle [Degrees]          Energy [eV]\n")
                    for i, j in zip(ref_x, data_delta):
                        fin.write("{:<20} {:<20}\n".format(i, j))

    # plot settings
    #marker = "."
    marker = ""
    lstyle = "-"
    plt.xlabel(xlabel)
    plt.ylabel("Energy / eV")
    plt.title(title)

    # ab initio data
    ref_x, ref_y = read_normed_output(ref_file)
    # shift x values by 2*pi for nicer plotting
    #ref_x = [i + 360 if (-180 < i < -150) else i for i in ref_x]

    # only for dihedrals
    #ref_x = [i+180 if i < 0 else i - 180 for i in ref_x]
    #ref_x = [180 + (-180 + i) if i < 0 else i for i in ref_x]
    #ref_x = [-180 - (180 - i) if i > 0 else i for i in ref_x]

    if x_offset is not None:
        plt.plot([i + x_offset if i < 0 else i for i in ref_x], ref_y, linestyle=lstyle, linewidth=0.2, marker=marker, label="ab initio")
    else:
        #plt.plot(ref_x, ref_y, linestyle="--", linewidth=0.02, marker=marker, label="ab initio")
        plt.plot(ref_x, ref_y, linestyle=lstyle, linewidth=1.0, marker=marker, label="ab initio")

    # plot data
    plot_data_files()
    #plt.legend(frameon=False)
    #plt.legend(edgecolor="white")
    #plt.legend(bbox_to_anchor=(0.5, 1.5))
    plt.legend(loc='upper right', bbox_to_anchor=(0.7, 1.0))
    plt.savefig("file.png")
    fig.show()
    #plt.show()


#def plot_results(rsgetters, param_dict, xlabel=r"Distance / $\AA$",
#                 ylabel=r"Energy / eV", title="Default", filename="default.png"):
#    """
#    Plot a list of Resultgetter instances
#    """
#    fig = plt.figure()
#    plt.xlabel(xlabel)
#    plt.ylabel(ylabel)
#    plt.title(title)
#
#    for rsgetter in rsgetters:
#        pdb.set_trace()
#        xvals, yvals = zip(*rsgetter.normed_results)
#        plt.plot(xvals, yvals, **param_dict)
#
#    plt.legend(loc='upper right', bbox_to_anchor=(0.7, 1.0))
#    plt.axhline(0, color='black', linestyle="--", linewidth=0.5)
#    fig.savefig(filename, dpi=600)
#    #fig.show()

#def plot_results(ax, rsgetters, param_dict):
#    """
#    Plot a list of Resultgetter instances
#    """
#    # horizontal line at 0.0
#    ax.axhline(0, color='black', linestyle="--", linewidth=0.5)
#
#    for rsgetter in rsgetters:
#        xvals, yvals = zip(*rsgetter.normed_results)
#        ax.plot(xvals, yvals, **param_dict)
#
#    ax.legend(loc='upper right', bbox_to_anchor=(0.7, 1.0))


def plot_multiple_histo(
        data_dict,
        xlim=(3.2, 6),
        xlabel=r"$\mathrm{r_{cog-cog}}$ / $\mathrm{\AA}$",
        png_name="default.png",
        x_ticks=(0.1, 0.5),
        y_ticks=(0.5, 1.0)):
    """
    Plot histograms for cog-cog distances.

    Parameters
    ----------
    **data : dict {str : list of floats}
        the data from each polymorph to plot the histogram from

    """
    # number of subplots depends on the number of values
    fontsize = 11
    nrows = len(data_dict)
    fig, axes = plt.subplots(nrows, 1, sharex=True, sharey=True, figsize=(10, 10))
    ax = fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    ax.set_xlabel(xlabel, fontweight='bold', fontsize=fontsize)
    ax.set_ylabel("Frequency (normalized)", fontweight='bold', fontsize=fontsize)

    for idx, (key, data) in enumerate(data_dict.items()):
        mu = np.mean(data)
        median = np.median(data)
        sigma = np.std(data)
        coeff_var = (sigma / mu) * 100
        num_bins = np.sqrt(len(data))
        num_bins = int(num_bins)

        if num_bins > 200:
            num_bins = 200  # cap number of bins to max. 200

        # make bins of same width
        data = sorted(data)
        x_min = data[0]
        x_max = data[-1]
        x_range = np.linspace(x_min, x_max, num_bins)

        # measured histogram
        nhist, bins, patches = axes[idx].hist(data, bins=x_range, density=True, histtype="step", align="mid", color=pantone_colors_hex["blue"])

        # fitted gaussian curve
        fitted_x, fitted_y = gnuplot_gaussfit(bins[1:], nhist)

        if idx == 0:
            clabel = "Fitted data"
        else:
            clabel = None

        axes[idx].plot(fitted_x, fitted_y, color=pantone_colors_hex["red"], linewidth=1.0, label=clabel)

        # line that goes through mu (=max. of gauss-function)
        #axes[idx].axvline(x=mu, linewidth=0.5, color="red", alpha=3.0, label="average")
        #axes[idx].axvline(x=median, linewidth=0.5, color="green", alpha=3.0, label="median")

        # coefficient of variation
        #coeff_var = (sigma / mu) * 100

        # Set a title for current subplot
        #axes[idx].set_title("Form {} $\mu={:.3f}$ $\sigma={:.3f}$ $VarCoeff={:.3f}$ %".format(key, mu, sigma, coeff_var))
        axes[idx].text(0.02, 0.97, "Form {} $\mu={:.3f}$ $\sigma={:.3f}$ $VarCoeff={:.3f}$ %".format(key, mu, sigma, coeff_var), transform=axes[idx].transAxes, fontsize=fontsize, fontweight='bold', )
        #axes[idx].text(0.02, 0.97, "Form {}".format(key), transform=axes[idx].transAxes, fontweight='bold')
        #axes[idx].legend(loc='upper right')

        # length of x-axis
        axes[idx].set_xlim(left=xlim[0], right=xlim[1])
        #axes[idx].set_ylim(0, 4.5)

        # ticker
        axes[idx].spines['right'].set_color('none')
        axes[idx].spines['top'].set_color('none')
        axes[idx].xaxis.set_minor_locator(ticker.MultipleLocator(x_ticks[0]))
        axes[idx].xaxis.set_major_locator(ticker.MultipleLocator(x_ticks[1]))
        axes[idx].yaxis.set_minor_locator(ticker.MultipleLocator(y_ticks[0]))
        axes[idx].yaxis.set_major_locator(ticker.MultipleLocator(y_ticks[1]))

    # distance between subplots
    fig.subplots_adjust(left=0.025, right=0.95, bottom=0.1, top=1.01)
    fig.tight_layout()
    fig.legend(frameon=False)
    plt.savefig(png_name, dpi=300)
    #plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ab_result_file",
                        help="Two column file with ab initio entities and values")
    parser.add_argument("md_result_files",
                        nargs="*",
                        help="Two column files with restrained md entities and values")
    parser.add_argument("-substract",
                        action="store_true",
                        help="Plot line which substracts ab initio reference from the shown ff iteration.")

    parser.add_argument("-x_offset",
                        type=float,
                        default=None,
                        help="Plot line which substracts ab initio reference from the shown ff iteration.")

    parser.add_argument("-xlabel",
                        default=None,
                        help="Title of x-axis")

    parser.add_argument("-title",
                        default=None,
                        help="Title of the plot")

    args = parser.parse_args()
    plot_all(args.ab_result_file, args.md_result_files,
             title=args.title, xlabel=args.xlabel,
             substract=args.substract, x_offset=args.x_offset)
