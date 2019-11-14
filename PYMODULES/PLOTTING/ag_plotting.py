
import pdb
import os
import collections
import pandas as pd
import statsmodels.api as sm
import scipy.stats
from scipy.optimize import curve_fit
import numpy as np
from ag_statistics import gnuplot_gaussfit


# check for xserver
import matplotlib as mpl

try:
    os.environ["DISPLAY"]
except KeyError:
    mpl.use('Agg')

import matplotlib.pyplot as plt


def gnuplot_gaussfit_plot(data, xlabel=None, output=None):
    """
    Plot the data as histogram and try fitting a gaussian curve over it.

    Parameters
    ----------
    data : list or numpy.ndarray
        Input values as single array.
    output : str
        Name of the output-file(s)

    """
    fig = plt.figure(figsize=(12, 8), facecolor='1.0')
    mu = np.mean(data)
    median = np.median(data)
    sigma = np.std(data)
    coeff_var = abs((sigma / mu) * 100)
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
    nhist, bins, _ = plt.hist(data, bins=x_range, density=True, histtype="step", align="mid")

    # fitted gaussian curve
    fitted_x, fitted_y = gnuplot_gaussfit(bins[1:], nhist)

    plt.plot(fitted_x, fitted_y, linewidth=1.0, label="Fitted data")
    plt.xlabel(xlabel, size=12)
    plt.ylabel("Frequency", size=12)
    plt.title(r"median={} $\mu$={:.3f} $\sigma$={:.3f} VarCoeff={:.3f} %".format(median, mu, sigma, coeff_var), size=16)
    fig.tight_layout()
    fig.legend(frameon=False)
    plt.savefig("{}.png".format(output), dpi=300)
