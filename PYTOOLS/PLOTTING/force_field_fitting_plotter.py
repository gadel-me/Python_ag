import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
import ag_statistics as ags
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


def plot_all(ref_file, data_files, xlabel=None, title=None):
    """
    Plot all data.

    Parameters
    ----------
    ref_file : str
    data_files : list of str
    xlabel : None or str

    """
    def plot_data_files():
        """
        Plot the given data and calculating chi square.

        Calculate chi square error and extrapolate the data to the reference.
        """
        # colormap
        norm = matplotlib.colors.Normalize(vmin=0, vmax=len(data_files))
        cmap = matplotlib.cm.get_cmap('Dark2')

        for index, data_file in enumerate(data_files):
            iteration = int(re.findall(r"\d+", data_file)[-1])
            data_x, data_y = read_normed_output(data_file)
            interp_y = np.interp(ref_x, data_x, data_y)
            chi_square_error = ags.chi_square_error(ref_y, interp_y)
            label = r"{:> 6} - $\chi^2$-err:{:> 8.7f} eV".format(iteration, chi_square_error)
            color = cmap(norm(index))
            plt.plot(data_x, data_y, linestyle="--", marker=marker, linewidth=0.5, markersize=0.5, color=color)
            plt.plot(ref_x, interp_y, linestyle="--", marker=marker, color=color, label=label)

            # testing
            if len(data_files) == 1:
                data_x = [i + 360 if (-180 < i < -150) else i for i in data_x]
                data_x2 = [i + 2.0 for i in data_x]
                plt.plot(data_x2, data_y, linestyle="--", marker=marker, linewidth=0.5, markersize=0.5, color="pink", label="shifted data")

    # plot settings
    marker = "."
    plt.figure()
    plt.xlabel(xlabel)
    plt.ylabel("Energy / eV")
    plt.title(title)

    # ab initio data
    ref_x, ref_y = read_normed_output(ref_file)
    # shift x values by 2*pi for nicer plotting
    ref_x = [i + 360 if (-180 < i < -150) else i for i in ref_x]
    plt.plot(ref_x, ref_y, linestyle="-", linewidth=2.2, marker=marker, label="ab initio")

    # plot data
    plot_data_files()
    plt.legend(bbox_to_anchor=(0.5, 0.5))
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ab_result_file",
                        help="Two column file with ab initio entities and values")
    parser.add_argument("md_result_files",
                        nargs="*",
                        help="Two column files with restrained md entities and values")
    parser.add_argument("-title",
                        default=None,
                        help="Title of the plot")

    args = parser.parse_args()
    plot_all(args.ab_result_file, args.md_result_files, title=args.title)
