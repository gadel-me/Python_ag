import re
import numpy as np
import itertools
import matplotlib
matplotlib.use('Agg')

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

import argparse
import ag_statistics as ags
from restrain_scan import norm_energy
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
    fig = plt.figure(figsize=(15, 10), dpi=300)

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
    marker = "."
    #marker = ""
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
    plt.show()


def plot_results(rsgetters, rslabels, linestyles=(), xlabel=r"Distance / $\AA$",
                 ylabel=r"Energy / eV", title="Default"):
    """
    Plot a list of Resultgetter instances
    """
    fig = plt.figure()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    for rsgetter, rslabel, linestyle in itertools.izip_longest(rsgetters, rslabels, linestyles):
        xvals, yvals = zip(*rsgetter.normed_results)

        if linestyle is None:
            plt.plot(xvals, yvals, label=rslabel, linestyle="-")
        else:
            plt.plot(xvals, yvals, label=rslabel, linestyle=linestyle)

    plt.legend(loc='upper right', bbox_to_anchor=(0.7, 1.0))
    plt.axhline(0, color='black', linestyle="--", linewidth=0.5)
    fig.show()


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
