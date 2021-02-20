import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.stats as stats
from scipy.stats import norm
import numpy as np
import math
import scipy.stats
from statsmodels.graphics.gofplots import qqplot
from statsmodels.graphics.tsaplots import plot_acf
import ag_statistics as ags
import pdb


def ask4frame(first_last, num_frames):
    """
    Ask for first or last frame.
    first_last      str; "first" or "last"
    """
    print(
        "Choose the {} frame of thermodynamic data (frames available: 0 to {})".format(
            first_last, num_frames
        )
    )
    return input("> ")


def ask4keyword(xory, keywords):
    """
    Ask for a keyword out of keywords.
    xory    str; "x-values" or "y-values"s
    """
    print(
        "Choose one thermo-keyword for the {} (keywords available: {}".format(
            xory, ", ".join(sorted(keywords))
        )
    )
    return input("> ")


def plot_qq(data, key):
    """
    Quantile-quantile plot.

    This plot generates its own sample of the idealized distribution that we are comparing with,
    in this case the Gaussian distribution. The idealized samples are divided into groups (e.g. 5),
    called quantiles. Each data point in the sample is paired with a similar member from the
    idealized distribution at the same cumulative distribution.

    The resulting points are plotted as a scatter plot with the idealized value on the x-axis and
    the data sample on the y-axis.

    A perfect match for the distribution will be shown by a line of dots on a 45-degree angle from
    the bottom left of the plot to the top right. Often a line is drawn on the plot to help make
    this expectation clear. Deviations by the dots from the line shows a deviation from the
    expected distribution.
    """
    data = np.array(data)
    # pdb.set_trace()
    qqplot(data, line="s")
    # Set a title for current subplot
    plt.title("QQ-Plot {}".format(key), fontweight="bold", fontsize=11)
    plt.show()


def plot_autocorrelation_function(data, key):
    """
    Plot the autocorrelation function.

    Sources: https://machinelearningmastery.com/gentle-introduction-autocorrelation-partial-autocorrelation/
    """
    plot_acf(data, title="Autocorellation Function {}".format(key))
    plt.show()


def plot_histogram(data, key, label=None):
    """
    Plot a histogram using given data
    """
    alpha = 0.05
    # data = stats.zscore(data)

    # Shapiro-Wilk Test (only for ~ 2000 samples)
    if len(data) <= 2000:
        print("Shapiro-Wilk Test")
        stat_shapiro, p_shapiro = stats.shapiro(data)
        normal_shapiro = p_shapiro > alpha

        if normal_shapiro:
            print("Sample looks Gaussian (fail to reject H0)\n")
        else:
            print("Sample does not look Gaussian (reject H0)\n")

    # D'Agostino's K^2 Test
    print("D'Agostino's K^2 Test")
    stat_agostino, p_agostino = stats.normaltest(data)
    normal_agostino = p_agostino > alpha

    if normal_agostino:
        print("Sample looks Gaussian (fail to reject H0)\n")
    else:
        print("Sample does not look Gaussian (reject H0)\n")

    # Anderson-Darling Test
    print("Anderson-Darling Test")
    result_anderson = stats.anderson(data, dist="norm")
    print("Statistic: {}".format(result_anderson.statistic))

    for idx in range(len(result_anderson.critical_values)):
        sl_anderson = result_anderson.significance_level[idx]
        cv_anderson = result_anderson.critical_values[idx]

        # check if the null hypothesis can be rejected (H0: normal distributed)
        normal_anderson = result_anderson.statistic < cv_anderson

        if normal_anderson:
            print(
                "{:> .3f}: {:> .3f}, data looks normal (fail to reject H0)".format(
                    sl_anderson, cv_anderson
                )
            )
        else:
            print(
                "{:> .3f}: {:> .3f}, data does not look normal (reject H0)".format(
                    sl_anderson, cv_anderson
                )
            )

    print("\n")

    # Compute the skewness of a data set, For normally distributed data,
    # the skewness should be about 0
    skewness = stats.skew(data)
    print("Skewness (should be 0): {:> 4.12f}".format(skewness))

    # determine if the skewness is close enough to 0 (statistically speaking)
    stat_skewness, p_skewness = stats.skewtest(data)
    print("Statistic (Skewness): {:> 4.12f}".format(stat_skewness))
    print("P-Value (skewness): {:> 4.12f}".format(p_skewness))
    print("\n")

    # Compute the kurtosis (tailing) of a data set, For normally distributed data,
    # the kurtosis should be about 0
    kurtosis = stats.kurtosis(data)
    print("Kurtosis (should be 0): {}".format(kurtosis))

    # determine if the kurtosis is close enough to 0 (statistically speaking)
    if len(data) > 20:
        stat_kurtosis, p_kurtosis = stats.kurtosistest(data)
        print("Statistic (Kurtosis): {:> 4.12f}".format(stat_kurtosis))
        print("P-Value (Kurtosis): {:> 4.12f}".format(p_kurtosis))
        print("\n")
    else:
        print("Need more than 20 values for Kurtosis-Test!")

    # /// define mu and sigma
    mu = np.mean(data)  # mean of distribution
    median = np.median(data)
    sigma = np.std(data)  # standard deviation of distribution

    # /// define bins
    num_bins = np.sqrt(len(data))
    num_bins = int(num_bins)  # only int

    if num_bins > 200:
        num_bins = 200  # cap number of bins to max. 200

    # make bins of same width
    data = sorted(data)
    x_min = data[0]
    x_max = data[-1]
    x_range = np.linspace(x_min, x_max, num_bins)

    plt.xlabel("delta", fontsize=10, fontweight="bold")
    plt.ylabel("Probability density", fontsize=10, fontweight="bold")

    # histogram
    nhist, bins, patches = plt.hist(
        data, bins=x_range, density=True, color="#E6E6FA", align="mid", cumulative=False
    )

    # define best fit line for given x_range, mu and sigma
    pdf_x_values = bins
    # pdf_y_values = mlab.normpdf(pdf_x_values, mu, sigma)
    pdf_y_values = norm.pdf(x_range)
    # pdb.set_trace()
    plt.plot(
        pdf_x_values,
        pdf_y_values,
        "g--",
        linewidth=0.5,
        alpha=0.5,
        label="Fitted data - estimated values",
    )

    fitted_x, fitted_y = ags.gnuplot_gaussfit(bins[1:], nhist)
    plt.plot(
        fitted_x,
        fitted_y,
        "r--",
        linewidth=0.5,
        alpha=0.5,
        label="Fitted data - gnuplot fitting",
    )

    # line that goes through mu (=max. of gauss-function)
    plt.axvline(x=mu, linewidth=0.5, color="red", alpha=3.0, label="average")
    plt.axvline(x=median, linewidth=0.5, color="green", alpha=3.0, label="median")

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)

    # coefficient of variation
    coeff_var = (sigma / mu) * 100

    # Set a title for current subplot
    plt.title(
        "Histogram of %r $\mu=%.4f$ $\sigma=%.4f$ $VarCoeff=%.4f$ %%"
        % (key, mu, sigma, coeff_var),
        fontweight="bold",
        fontsize=11,
    )

    # chi square test to check the quality of the fit
    histo, bin_edges = np.histogram(data, bins=num_bins, normed=False)
    a1, b1 = stats.norm.fit(data)
    cdf = stats.norm.cdf(bin_edges, a1, b1)
    # scaling_factor = len(data) * (x_max - x_min) / num_bins
    scaling_factor = len(data)
    # expected frequencies (haufigkeiten)
    expected_values = scaling_factor * np.diff(cdf)
    chisquare_results = stats.chisquare(histo, f_exp=expected_values, ddof=2)
    print(chisquare_results, "\n")
    # pdb.set_trace()

    # Kolmogorov-Smirnov test for goodness of fit
    kstest_results = stats.kstest(data, "norm")
    print(kstest_results, "\n")

    # z_scores = stats.zscore(data)
    # Kolmogorov-Smirnov test for goodness of fit
    # zscore_kstest_results = stats.kstest(z_scores, "norm")
    # print(zscore_kstest_results, "\n")

    # if kstest_results[1] > 0.05:
    #    print("{} > 0.05: Normal distribution seems identical to given distribution (failed to reject H0)\n".format(kstest_results[1]))
    # else:
    #    print("{} < 0.05: Normal distribution is not identical to given distribution (reject H0)\n".format(kstest_results[1]))

    plt.legend()
    plt.show()


def plot_xy(xvals, yvals, xkey, ykey, linreg=False):
    """
    Plot y-values against x-values.

    Plot the x against the y values.
    """
    # window layout#
    xyfig = plt.figure()
    xyfig.canvas.set_window_title("Scatter-Plot %s vs. %s" % (xkey, ykey))
    plt.xlabel(xkey, fontsize=12, fontweight="bold")
    plt.ylabel(ykey, fontsize=12, fontweight="bold")
    plt.legend(
        bbox_to_anchor=(0.0, 1.0, 1.0, 0.1),
        loc="upper right",
        borderaxespad=0.0,
        frameon=True,
        shadow=True,
        numpoints=1,
        prop={"size": 10},
    )
    plt.title("Scatter-Plot: %r vs. %r" % (xkey, ykey), fontweight="bold")

    if linreg is True:
        slope, intercept, r_value, p_value, std_err = stats.linregress(xvals, yvals)
        print("{:<10s} {:> 4.12f}".format("Slope:", slope))
        print("{:<10s} {:> 4.12f}".format("Intercept:", intercept))
        print("{:<10s} {:> 4.12f}".format("R:", r_value))
        print("{:<10s} {:> 4.12f}".format("P:", p_value))
        print("{:<10s} {:> 4.12f}".format("Std.err.:", std_err))

    plt.plot(xvals, yvals, "r-", antialiased=True)
    plt.show()
