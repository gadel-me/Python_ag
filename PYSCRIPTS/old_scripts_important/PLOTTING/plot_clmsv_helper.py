from __future__ import print_function, division
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats as stats
import numpy as np
import math


def ask4frame(first_last, num_frames):
    """
    Ask for first or last frame.
    first_last      str; "first" or "last"
    """
    print("Choose the {} frame of thermodynamic data (frames available: 0 to {})".format(first_last, num_frames))
    return raw_input("> ")


def ask4keyword(xory, keywords):
    """
    Ask for a keyword out of keywords.
    xory    str; "x-values" or "y-values"s
    """
    print("Choose one thermo-keyword for the {} (keywords available: {}".format(xory, ", ".join(sorted(keywords))))
    return raw_input("> ")


def plot_histogram(data, key, label=None):
    """
    Plot a histogram using given data
    """
    # /// define mu and sigma
    mu    = np.mean(data)  # mean of distribution
    sigma = np.std(data)   # standard deviation of distribution

    # /// define bins
    num_bins = np.sqrt(len(data))
    num_bins = math.ceil(num_bins)  # only int

    if num_bins > 200:
        num_bins = 200  # cap number of bins to max. 200

    # make bins of same width
    data = sorted(data)
    x_min = data[0]
    x_max = data[-1]
    x_range = np.linspace(x_min, x_max, num_bins)
    plt.xlabel("delta", fontsize=10, fontweight='bold')
    plt.ylabel("Probability density", fontsize=10, fontweight='bold')

    # histogram
    n, bins, patches = plt.hist(data, bins=x_range, normed=True, color="#E6E6FA",
                                align="mid", cumulative=False)

    # define best fit line for given x_range, mu and sigma
    pdf_x_values = bins
    pdf_y_values = mlab.normpdf(pdf_x_values, mu, sigma)
    plt.plot(pdf_x_values, pdf_y_values, "r--", linewidth=0.5)

    # line that goes through mu (=max. of gauss-function)
    plt.axvline(x=mu, linewidth=0.5, color="red", alpha=1.0)

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)

    # Set a title for current subplot
    plt.title("Histogram of %r $\mu=%.4f$ $\sigma=%.4f$" % (key, mu, sigma),
              fontweight='bold', fontsize=11)
    plt.show()


def plot_xy(xvals, yvals, xkey, ykey, linreg=False):
    """
    Plot y-values against x-values
    """
    # window layout#
    xyfig = plt.figure()
    xyfig.canvas.set_window_title("Scatter-Plot %s vs. %s" % (xkey, ykey))
    plt.xlabel(xkey, fontsize=12, fontweight='bold')
    plt.ylabel(ykey, fontsize=12, fontweight='bold')
    plt.legend(bbox_to_anchor=(0., 1., 1., .1), loc="upper right",
               borderaxespad=0., frameon=True, shadow=True, numpoints=1,
               prop={'size': 10})
    plt.title("Scatter-Plot: %r vs. %r" % (xkey, ykey), fontweight='bold')

    if linreg is True:
        slope, intercept, r_value, p_value, std_err = stats.linregress(xvals, yvals)

    plt.plot(xvals, yvals, "r-", antialiased=True)
    plt.show()
