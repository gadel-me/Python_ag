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
import scipy.stats
from scipy.optimize import curve_fit
import numpy as np
#import pdb


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


def test_gauss_shape(test, data, min_val=-0.3, max_val=0.3, filename=None):
    """
    Carry out the Kurtosis, Skewness test or a QQ-Plot.

    Tests the shape of the gaussian curve. A value of zero means the curve
    has a total gaussian shape. Values beyond min and max indicate that
    the curve is quite far away from a symmetric gaussian shape.

    Sources: https://en.wikipedia.org/wiki/Kurtosis
             https://en.wikipedia.org/wiki/Skewness

    Input:
        > test      str;
                    skewness/kurtosis/qq
        > min_val   float;
                    minimum, skewness/kurtosis should not be smaller than this
        > max_val   float;
                    maximum, skewness/kurtosis should not be larger than this

    Returns:
        bool; True (skewness or kurtosis is between min and max) or False (is not)

    Raises:
        Warning; less than 20 values given for kurtosis
        NameError; wrong name for test

    """
    if test == "skewness" or test == "kurtosis":
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

        result_str = "minimum: {1}; {0}: {2: .3f}; maximum: {3}".format(test, min_val, value,
                                                                        max_val)
        result = min_val <= value <= max_val
    elif test == "qq":
        value = scipy.stats.probplot(data, dist="norm", fit=True)
        result_str = "QQ - Straight line\nslope: {}, interecept: {}, r^2: {}".format(value[1][0], value[1][1],
                                                                                     value[1][2]**2)
        result = value[1][2]**2 > 0.99
    else:
        raise NameError("Test named '{}' not known!".format(test))

    # print results to file or console
    if filename is not None:
        with open(filename, "a") as f_open:
            f_open.write(result_str)
    else:
        print(result_str)

    return result


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

    ### Gaussian function to fit with
    def gauss(x, *p):
        A, mu, sigma = p
        return A / np.sqrt(2 * sigma**2 * np.pi) * np.exp(-(x - mu)**2 / (2.0 * sigma**2))

    ### RMSE to evaluate goodness of fit
    def rmse(predictions, targets):
        return np.sqrt(((predictions - targets) ** 2).mean())

    xmax = max(x_values)  # minimum x value of data to plot
    xmin = min(x_values)  # maximum x value of data to plot
    p0 = [1., (xmax + xmin) / 2, (xmax - xmin) / 8]  # initial value to start fitting with

    ### fitting procedure
    coeff, var_matrix = curve_fit(gauss, x_values, y_values, p0=p0, maxfev=1000000)

    hist_fit = gauss(x_values, *coeff)

    gauss_data = np.array((x_values, hist_fit))
    rmse_val = rmse(gauss_data[1], y_values)

    if debug is True:
        print("Rms error: {0}".format(rmse_val))
        print('Fitted mean = {0}'.format(coeff[1]))
        print('Fitted standard deviation = {0}\n'.format(coeff[2]))

    return(gauss_data)
