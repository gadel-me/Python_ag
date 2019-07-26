"""
Module to carry out statstical tests on a certain set of data.

All normality tests test against normal distribution, therefor they can only fail
to reject normality but they never accept normality.

Keep in mind  when dealing with big amounts of data even small deviations from
normality may lead to a false rejection of H0.

Source: https://machinelearningmastery.com/a-gentle-introduction-to-normality-tests-in-python/
        https://www.rdocumentation.org/packages/TeachingDemos/versions/2.10/topics/SnowsPenultimateNormalityTest
        https://stackoverflow.com/questions/7781798/seeing-if-data-is-normally-distributed-in-r/7788452#7788452
        https://medium.com/@rrfd/testing-for-normality-applications-with-python-6bf06ed646a9
        https://www.crashkurs-statistik.de/vorgehen-bei-hypothesentests/
        http://onlinestatbook.com/2/advanced_graphs/q-q_plots.html
        http://desktop.arcgis.com/en/arcmap/latest/extensions/geostatistical-analyst/normal-qq-plot-and-general-qq-plot.htm
"""

from __future__ import print_function, division
import pdb
import os
import collections
import pandas as pd
import statsmodels.api as sm
import scipy.stats
from scipy.optimize import curve_fit
import numpy as np


# use matplotlib w / w/o xserver
import matplotlib as mpl

try:
    os.environ["DISPLAY"]
except KeyError:
    mpl.use('Agg')

import matplotlib.pyplot as plt


def _print_result(p_value, alpha, statistic, test_name, filename=None):
    """
    Print test results using this helper function.

    Parameters:
        > test_result   bool;
                        result of the normality test
        > p_value       float;
                        p-value to fail to reject H0
        > alpha         float;
                        threshold to reject/fail to reject H0
        > statistic     float;
                        statistical value of test
        > test_name     str;
                        name of the test that was carried out
        > filename      str;
                        optional; print or write results to a file

    """
    # strings to print
    fail_reject_h0 = "{: >.5f} (P-value) > {:>.5f}; Sample looks Gaussian (fail to reject H0)\n".format(
        p_value, alpha
    )
    reject_h0 = "{: >.5f} (P-value) < {:>.5f}; Sample does not look Gaussian (reject H0)\n".format(
        p_value, alpha
    )
    statistic_str = "Statistic: {:> .8f}".format(statistic)

    if filename is None:
        print(test_name)
        print(statistic_str)
    else:
        with open(filename, "a") as f_open:
            f_open.write(test_name)
            f_open.write("\n")
            f_open.write(statistic_str)
            f_open.write("\n")

    if p_value >= alpha:

        if filename is None:
            print(fail_reject_h0)
        else:
            with open(filename, "a") as f_open:
                f_open.write(fail_reject_h0)

    else:

        if filename is None:
            print(reject_h0)
        else:
            with open(filename, "a") as f_open:
                f_open.write(reject_h0)


def test_normality(test, data, alpha=0.05, filename=None):
    """
    Carry out a Shapiro-Wilk, Agostino K^2, Anderson-Darling, ... Test.

    Check if the given data is not normally distributed. The null hypothesis and
    the alternative hypothesis are defined as following:

        > H0: data is normally distributed
        > H1: data is not normally distributed


    To interpret the p-value a threshold level called alpha is chosen which is
    typically 5% (or 0.05). It applies as following:
        > p-value > 0.05: Failed to reject H0
        > p-value < 0.05: Reject H0

    The Shapiro-Wilk Test tests if H0 can or cannot be rejected and if H1 has
    to be accepted and only works for a sample size smaller than 2000.

    Input:
        > test          str; name of the test to be carried out. The following
                        are allowed: "shapiro", "agostino", "anderson",
                        "skewness", "kurtosis"
        > data          list or np-array; the data set to test
        > alpha         float; probability threshold
        > filename      None or str; if a file name is given, write results to
                        a file, else print the results to the screen

    Returns:
        bool; True (cannot reject H0) or False (reject H0)

    Raises:
        > Warning; shapiro with more than 2000 samples
        > NameError; wrong name for test given

    """
    if test == "shapiro":
        test_name = "Shapiro-Wilk"

        # abort if > 2000 values
        if len(data) > 5000:
            raise Warning("More than 2000 values given!")

        stat, p_value = scipy.stats.shapiro(data)
    elif test == "agostino":
        test_name = "D'Agostino's K^2"
        stat, p_value = scipy.stats.normaltest(data)
    elif test == "anderson":
        # returns true or false since it has several levels of significance
        test_name = "Anderson-Darling"
        result_anderson = scipy.stats.anderson(data, dist="norm")
        stat = result_anderson.statistic

        for idx in range(len(result_anderson.critical_values)):
            alpha_anderson = result_anderson.significance_level[idx]
            cv_anderson = result_anderson.critical_values[idx]

            # check if the null hypothesis can be rejected (H0: normal distributed)
            #_print_result(cv_anderson, stat, stat, test_name, filename)

            # if H0 is ok at any confidence level, stop further testing
            if (cv_anderson > result_anderson.statistic) and (alpha_anderson == alpha):
                return cv_anderson > result_anderson.statistic

        # anderson darling is false if it rejected H0 for all significance levels
        return False

    elif test == "skewness":
        test_name = "Skewness"
        stat, p_value = scipy.stats.skewtest(data)
    elif test == "kurtosis":
        test_name = "Kurtosis"
        stat, p_value = scipy.stats.kurtosistest(data)
    elif test == "chi_square":
        pass
    else:
        raise NameError("Test named '{}' not known!".format(test))

    _print_result(p_value, alpha, stat, test_name, filename)
    return p_value > alpha


def _write_data(data1, data2, output):
    """
    Write the data into two columns of a csv file.

    This is just a helper function for qq_test and probplot_test.
    """
    summary = collections.OrderedDict({"Experimental Quantiles": data1, "Theoretical Quantiles": data2})
    table = pd.DataFrame(summary)
    table.to_csv("{}.csv".format(output), index=False)
    return table


def _write_plot(data1, data2, output, slope, intercept, title):
    """
    Save the plot of data1 and data2.

    This is just a helper function for qq_test and probplot_test.
    """
    qq_fig = plt.figure(figsize=(12, 8), facecolor='1.0')
    plt.plot(data1, data2, "o")
    plt.plot(data1, intercept + slope * data1, "r--", label="fitted line")
    plt.title(title, size=24)
    plt.xlabel("Theoretical quantiles", size=12)
    plt.ylabel("Experimental quantiles", size=12)
    plt.tick_params(labelsize=10)
    plt.savefig("{}.png".format(output))
    #plt.show()
    plt.close(qq_fig)


def qq_test(data, rsquare_thrsh=0.99, output=None, save_plot=False):
    """
    Carry out a Quantile-Quantile plot with the experimental data and the theoretical gaussian shaped data.

    "Many standard statistical procedures require normally distributed data.
    One way to assess if your data is normally distributed is quantile-quantile
    plot or q-q plot. In this approach quantiles of a tested distribution are
    plotted against quantiles of a known distribution as a scatter plot.
    If distributions are similar the plot will be close to a straight line.
    We will plot our data against a normal distribution to test if our data
    is distributed normally." (1)
    Returns the r^2 value from the linear regression. If r^2 > 'rsquare_thrsh' then
    the data should be equal to a normal distribution.

    Sources
    -------
    (1) https://medium.com/@rrfd/testing-for-normality-applications-with-python-6bf06ed646a9
    (2) https://scientificpythonsnippets.com/index.php/2-uncategorised/6-q-q-plot-in-python-to-test-if-data-is-normally-distributed

    Parameters
    ----------
    data : list
        the data to check

    rsquare_thrsh : float or int, default: 0.99
        rsquare of the linear regression must be at least

    output : str or None, optional
        output name of the files that will be created

    save_plot : bool, optional
        Save a picture of the plot if True

    Returns
    -------
    rsquare >= rsquare_thrsh : bool

    """
    sorted_data = sorted(data)
    # created gaussian shaped data 'norm' from given data
    norm = np.random.normal(0, 2, len(sorted_data))
    # for Q-Q plot all data must be sorted
    sorted_norm = np.sort(norm)
    # trend line from linear regression
    slope, intercept, rvalue, pvalue, stderr = scipy.stats.linregress(sorted_norm, sorted_data)
    rsquare = rvalue**2

    # write raw data to file
    if output is not None:
        _write_data(sorted_data, sorted_norm, output)

    # save the plot as png file
    if save_plot is True:
        _write_plot(sorted_norm, sorted_data, output, slope, intercept, "Q-Q plot - " + r"$r^2$" + "={:> 2.4f}".format(rsquare))

    return rsquare >= rsquare_thrsh


def probability_plot(data, rsquare_thrsh=0.99, output=None, save_plot=False):
    """
    Carry out a probability plot (not Q-Q and not P-P) with the experimental data.

    "probplot generates a probability plot, which should not be confused with
    a Q-Q or a P-P plot. Statsmodels has more extensive functionality of this
    type, see statsmodels.api.ProbPlot." (1)

    Sources
    -------
    (1) https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.probplot.html

    Parameters
    ----------
    data : list
        the data to check

    rsquare_thrsh : float or int, default: 0.99
        rsquare of the linear regression must be at least

    output : str or None, optional
        output name of the files that will be created

    save_plot : bool, optional
        Save a picture of the plot if True

    Returns
    -------
    rsquare >= rsquare_thrsh : bool

    """
    (osm, osr), (slope, intercept, rvalue) = scipy.stats.probplot(data, dist="norm", fit=True)
    rsquare = rvalue**2

    # write raw data to file
    if output is not None:
        _write_data(osr, osm, output)

    # save the plot as png file
    if save_plot is True:
        _write_plot(osm, osr, output, slope, intercept, r"Probability plot ($\r^2$={})".format(rsquare))

    return rsquare >= rsquare_thrsh


def test_gauss_shape(test, data, min_val=-0.3, max_val=0.3):
    """
    Carry out the Kurtosis and Skewness test.

    Tests the shape of the gaussian curve. A value of zero means the curve
    has a total gaussian shape. Values beyond min and max indicate that
    the curve is quite far away from a symmetric gaussian shape.

    Parameters
    ----------
    min_val : float
        Minimum; skewness/kurtosis should not be smaller than this
    max_val : float
        Maximum; skewness/kurtosis should not be larger than this

    Returns
    --------
    min_val <= value <= max_val
        True (skewness or kurtosis is between min and max) or False (is not)

    Raises
    ------
    Warning; less than 20 values given for kurtosis
    NameError; wrong name for test

    Sources
    -------
    (1) https://en.wikipedia.org/wiki/Kurtosis
    (2) https://en.wikipedia.org/wiki/Skewness

    """
    # kurtosis only valid if sample is larger than 20
    if test == "skewness":
        value = scipy.stats.skew(data)
    elif test == "kurtosis":

        if len(data) > 20:
            value = scipy.stats.kurtosis(data)
        else:
            raise Warning("Kurtosis needs at least 20 values!")

    else:
        raise NameError("Test named '{}' not known!".format(test))

    result_str = "minimum: {1}; {0}: {2: .3f}; maximum: {3}".format(test, min_val, value, max_val)
    print(result_str)
    return min_val <= value <= max_val


def chi_square_error(data1, data2):
    """
    Calculate the chi square error of two sets of data.

    Parameters
    ----------
    data1 : list/np-array of floats
        first set of data to compare

    data2 : list/np-array of floats
        second set of data to compare

    Returns
    -------
    sum(chi_squares) : float
        the chi square value

    Sources: http://cmt.dur.ac.uk/sjc/thesis_dlc/node88.html
    """
    chi_squares = []

    for en1, en2 in zip(data1, data2):
        chi_squares.append(abs(en1 - en2)**2)

    return sum(chi_squares)


def gnuplot_gaussfit(x_values, y_values, debug=False):
    """
    Fit a gaussian curve for the data given.

    Parameters
    ----------
    x_values : array_like
        1-d array with data

    y_values : array_like
        1-d array with data
    """
    # Gaussian function to fit with
    def gauss(x, *p):
        """
        Equation of the gaussian function.
        """
        A, mu, sigma = p
        return A / np.sqrt(2 * sigma**2 * np.pi) * np.exp(-(x - mu)**2 / (2.0 * sigma**2))

    # RMSE to evaluate goodness of fit
    def rmse(predictions, targets):
        return np.sqrt(((predictions - targets) ** 2).mean())

    xmax = max(x_values)  # minimum x value of data to plot
    xmin = min(x_values)  # maximum x value of data to plot
    p0 = [1., (xmax + xmin) / 2, (xmax - xmin) / 8]  # initial value to start fitting with

    # fitting procedure
    coeff, var_matrix = curve_fit(gauss, x_values, y_values, p0=p0, maxfev=1000000)
    hist_fit = gauss(x_values, *coeff)
    gauss_data = np.array((x_values, hist_fit))
    rmse_val = rmse(gauss_data[1], y_values)

    if debug is True:
        print("Rms error: {0}".format(rmse_val))
        print('Fitted mean = {0}'.format(coeff[1]))
        print('Fitted standard deviation = {0}\n'.format(coeff[2]))

    return(gauss_data)
