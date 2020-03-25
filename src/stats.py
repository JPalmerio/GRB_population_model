import time
import logging
import numpy as np
from scipy.stats import mstats, ks_2samp
from astroML.utils import check_random_state
from plotting_functions import plot_CDFs_and_KS_results
from io_grb_pop import read_column

log = logging.getLogger(__name__)


def draw_from_cdf_file(filename, N_draws, **args):
    """
        Draw from an ascii file that contains two columns:
        x, CDF(x)
    """
    value_range = read_column(filename, 0)
    cdf = read_column(filename, 1)
    draws = np.random.rand(N_draws)
    values = value_range[cdf.searchsorted(draws)]
    return values


def create_pdf_from_cdf_file(filename, **args):
    x = read_column(filename, 0, **args)
    cdf = read_column(filename, 1, **args)
    pdf = np.zeros(len(cdf))
    pdf = cdf[1:]-cdf[:-1]
    for i in range(1,len(cdf)):
        pdf[i] = cdf[i]-cdf[i-1]
    return x, pdf


def pBIL(mod, obs, epsilon=1e-3, sum_ln_oi_factorial=None, verbose=False):
    """
        The parametric Bayesian Indirect Likelihood (pBIL) is a goodness
        of fit estimator based on a the idea that the number of objects
        in each bin follows a Poissonian law. The likelihood is then:
        likelihood = prod( exp(-s_i) * s_i**o_i / o_i! )
        ln(likelihood) = sum( o_i * ln(s_i) - s_i - ln(o_i!) )
        where s_i is the predicted number of objects in bin_i from the
        model and o_i is the observed number of objects in bin_i
        Epsilon is a user-defined constant added to empty bins to avoid
        problems when computing the log. It should be small and does
        not have an impact on the value of the likelihood as long as the
        model fits the observations fairly well (i.e. no empty bins)
    """
    # Replace empty bins
    mod_no_empty_bins = np.where(mod == 0, epsilon, mod)
    if sum_ln_oi_factorial is None:
        # Calculate sum( ln(o_i!) ) factor
        sum_ln_oi_factorial = 0.
        for i, val in enumerate(obs):
            if val != 0:
                sum_ln_oi_factorial += val*np.log(val) - val
        if verbose:
            log.info('In pBIL, sum_ln_oi_factorial = {}'.format(sum_ln_oi_factorial))
    return np.sum(obs * np.log(mod_no_empty_bins) - mod_no_empty_bins) - sum_ln_oi_factorial


def chi2(mod, obs, err):

    return np.sum((np.abs(mod-obs)/err)**2)


def unbinned_empirical_cdf(data, weights=1):
    """
        From the answer of Dave at http://stackoverflow.com/questions/3209362/how-to-plot-empirical-cdf-in-matplotlib-in-python
        Note : if you wish to plot, use arg drawstyle='steps-post', I found it is the most accurate
        description of the data

        Parameters:
        ------------

        data : [array]
            The data to convert into a Cumulative Distribution Function.

        weights : [array]
            The weights for the data.

        Returns:
        ---------

        sorted_data : [array]
            The array of the data sorted.

        CDF : [array]
            The cumulative distribution function that follows the formal definition of CDF(x) = "number of samples <= x"/"number of samples"
    """
    sorted_data = np.sort(data)
    if isinstance(weights, np.ndarray):
        # create 2D array with data and weights
        arr = np.column_stack((data, weights))
        # Sort them by ascending data, need to use this method and not np.sort()
        arr = arr[arr[:,0].argsort()]
        CDF = np.cumsum(arr[:,1]).astype(float)
        CDF /= CDF[-1]
    else:
        CDF = np.arange(1, len(sorted_data)+1) / float(len(sorted_data))
    return sorted_data, CDF


def subsample_and_KS(df1, df2, N_sub, key, confidence=95.0, N_bs=100, precision=500, bins=None, create_CDF=True,
    show_plot=False, label1=None, label2=None, subsample1=True, subsample2=True, color1=None, color2=None):
    """
        Compute the K-S test between a subsample of size N_sub.
        Repeat the K-S test for N_bs bootstraps and create a distribution of p-values.
    """

    t1 = time.time()
    if subsample1:
        subsamp1 = subsample_with_bootstrap(df1, key=key, N_sub=N_sub, N_bs=N_bs)
    else:
        subsamp1 = df1[key].to_numpy()[:,np.newaxis] * np.ones(N_bs)
        subsamp1 = subsamp1.T
    if subsample2:
        subsamp2 = subsample_with_bootstrap(df2, key=key, N_sub=N_sub, N_bs=N_bs)
    else:
        subsamp2 = df2[key].to_numpy()[:,np.newaxis] * np.ones(N_bs)
        subsamp2 = subsamp2.T
    t2 = time.time()
    log.info(f"Subsampling done in {t2-t1:.3f} s")

    t1 = time.time()
    D_stat, p_value = compute_KS_with_bootstrap(subsamp1, subsamp2, N_bs)
    t2 = time.time()
    log.info(f"KS calculations done in {t2-t1:.3f} s")

    if create_CDF:
        t1 = time.time()
        CDF1 = CDF_with_bootstrap(subsamp1, bins=bins)
        CDF2 = CDF_with_bootstrap(subsamp2, bins=bins)
        med1, lw1, up1 = compute_CDF_quantiles(CDF1, confidence=confidence)
        med2, lw2, up2 = compute_CDF_quantiles(CDF2, confidence=confidence)
        t2 = time.time()
        log.info(f"CDF calculation done in {t2-t1:.3f} s")

    pfrac = len(p_value[np.where(p_value < (1.-confidence/100.))[0]])/len(p_value)

    if show_plot:
        plot_CDFs_and_KS_results(bins=bins, med1=med1, med2=med2, lw1=lw1, lw2=lw2, up1=up1, up2=up2,
                                 D_stat=D_stat,
                                 p_value=p_value,
                                 confidence=confidence,
                                 pfrac=pfrac,
                                 label1=label1,
                                 label2=label2,
                                 color1=color1,
                                 color2=color2)

    return pfrac


def subsample_with_bootstrap(df, key, N_sub, N_bs=100):
    """
        Creates an 2D array of subsamples from df[key] of size (N_bs, N_sub)
        where N_bs in the number of bootstraps and N_sub is the size of the subsample
    """
    subsamples = np.zeros((N_bs, N_sub))
    for i in range(N_bs):
        # Create a subsample of size N_sub from the column 'key' of df and turn it into a numpy array
        subsamples[i] = df[key].sample(N_sub).to_numpy()
    return subsamples


def CDF_with_bootstrap(subsamples, bins):
    """
        Takes in a subsample array of size (N_bs, N_sub) and computes the CDF for each bootstrap realization
        The CDFs are sampled given a certain precision determined from the bins.
        N_bs stands for number of bootstraps
        N_sub stands for size of the subsample
    """
    N_bs = subsamples.shape[0]
    precision = bins.shape[0] - 1
    # Compute the CDF for each subsample realization
    CDF_real = np.zeros((N_bs, precision))
    for i in range(N_bs):
        hist, bins_ = np.histogram(subsamples[i], bins=bins)
        CDF_real[i,:] = np.cumsum(hist).astype(float)/float(np.sum(hist))

    return CDF_real


def compute_CDF_quantiles(CDFs, confidence=95.0):
    """
        Takes a 2D array of CDFs of size (N_bs, N_bins).
        N_bs stands for the number of bootstraps
        N_bins stands for the number of bins within the CDF.
        Returns the median, lower and upper bounds at the desired confidence level.
    """

    N_bins = CDFs.shape[1]
    # Create percentiles:
    lower_percentile = (1. - confidence/100.)/2.
    upper_percentile = 1. - lower_percentile
    # Compute the percentiles for each bin
    lower = np.zeros(N_bins)
    median = np.zeros(N_bins)
    upper = np.zeros(N_bins)
    for i in range(N_bins):
        q = mstats.mquantiles(CDFs[:,i], prob=[lower_percentile, 0.5, upper_percentile])
        median[i] = q[1]
        upper[i] = q[2]
        lower[i] = q[0]
    return median, lower, upper


def compute_KS_with_bootstrap(samples1, samples2, N_bs):
    """
        Compute the 2 sample K-S test N_bs times where
        N_bs stands for the number of bootstrap samples
        Return the D-statistic and the pvalue for each bootstrapped
        realization in the form of 1D arrays of size N_bs
    """
    D_stat = np.zeros(N_bs)
    p_value = np.zeros(N_bs)
    for i in range(N_bs):
        D_stat[i], p_value[i] = ks_2samp(samples1[i], samples2[i])
    return D_stat, p_value


def compute_CDF_bounds_by_MC(sample, sample_errp, sample_errm=None, sample_ll=None, sample_ul=None,
    ll_max_val=None, ul_min_val=None, weights=None, weights_err=None, weqs=False, confidence=95.0,
    bins=None, positive=False, precision=1000, precision_pdf=1000, N_MC=1000, bootstrap=False,
    verbose=False, random_state=None, show_plot=False, ax=None, color='k', show_median=True,
    **kwargs):
    """
        Function to compute the lower and upper bounds of a cumulative distribution function for quantities with errors.
        Note if this function takes too long, try reducing N_MC or precision.
        Returns the middle point for the bins, the median CDF computed, the lower and upper CDF envelopes and the figure on which it has plotted.

        Parameters:
        -----------

        sample : [numpy array]
            Numpy array of the data for which to compute the bounds.

        sample_errp : [numpy array]
            Numpy array of the positive error on the data for which to compute the bounds.

        sample_errm : [numpy array]
            Default is None
            Numpy array of the negative error on the data for which to compute the bounds.
            If None, the errors are assumed to be symetric and equal to sample_errp

        sample_ll : [numpy array]
            Default is None
            Numpy array of the lower limits on the data for which to compute the bounds.
            Expects 0 if not a limit, 1 if a limit.

        sample_ul : [numpy array]
            Default is None
            Numpy array of the upper limits on the data for which to compute the bounds.
            Expects 0 if not a limit, 1 if a limit.

        weights : [numpy array]
            Default is None
            Numpy array of the weights on the data for which to compute the bounds.
            If None, the weights are assumed to be equal and no weighting is computed.

        weights_err : [numpy array or list of 2 numpy arrays]
            Default is None
            Numpy array of the error on the weights.
            If a numpy array is provided, the weights are drawn from an gaussian with sigma equal to weights.
            If a list of 2 numpy arrays is provided, the first numpy array is assumed to be the plus error and the second the minus.
            ex : [weights_errp, weights_errm]

        weqs : [bool]
            Default is False.
            If True, the weights are assumed equal to the sample and are not redrawn?

        confidence : [float]
            Default is 95.0
            Confidence level of the bounds (in percent).

        bins : [numpy array]
            Default is None
            Bins in which to compute the cumulative distribution function.
            If None, will use the conservative assumption that the smallest bin value is the min value of the sample minus 5 time the maximum error of the sample (opposite for biggest bin value).

        positive : [bool]
            Default is False
            Set to True if the quantity can not be negative (a distance for example), will reject negative drawings.

        bootstrap : [bool]
            Default is False
            Set to True to perform bootstrap resampling.

        ll_max_val : [numpy array]
            Default is None
            Numpy array of the minimum value to use as a lower bound for the prior on lower limits.
            If lower limits are provided (with sample_ll) but ll_max_val is not specified, by default will use the maximum value of the sample plus 5 times the largest error

        ul_min_val : [numpy array]
            Default is None
            Numpy array of the maximum value to use as an upper bound for the prior on upper limits.
            If upper limits are provided (with sample_ul) but ul_min_val is not specified, by default will use the minimum value of the sample minus 5 times the largest error

        precision : [int]
            Default is 1000
            Numbers of bins in which to compute the cumulative distribution function.

        precision_pdf : [int]
            Default is 1000
            Numbers of points on which to sample the probability distribution function of each data point.

        N_MC : [int]
            Default is 1000
            Numbers of Monte Carlo realizations of the sample.

        show_plot : [bool]
            Default is False
            If True, will create plots.

        ax : [matplotlib axes instance]
            Default is None
            If specified, will plot the median and upper and lower bounds on the ax.

        color : [str]
            Default is 'k'
            Color of the CDF to plot.

        **kwargs : [dict]
            Any matplotlib key words to pass to the median plot.

        Returns
        -------
        bins_mid : [array]
            Numpy array of values for the middle of the bins used in the cumulative distribution function

        median : [array]
            Numpy array of values for the median at each bin.

        lower : [array]
            Numpy array of values for the lower bound at each bin.

        upper : [array]
            Numpy array of values for the upper bound at each bin.

        fig : [matplotlib Figure]
            The figure on which is drawn the plot. Returns None if no show_plot is False.

    """
    if show_plot:
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    # Check inputs
    if not isinstance(sample, np.ndarray):
        raise TypeError('sample must be numpy array.')

    if verbose:
        log.info("In compute_CDF_bounds_by_MC: initializing...")
    # If plot is demanded but no ax is given, create figure and plot also indivudal points' PDFs
    if show_plot and (ax is None):
        fig = plt.figure(figsize=(8,5))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)
        fig.subplots_adjust(hspace=0)
        plt.setp(ax1.get_xticklabels(),visible=False)
        ax1.tick_params(axis='x',which='both',bottom=False)
    elif show_plot:
        ax1 = ax
        fig = plt.gcf()
    else:
        fig = None

    sample_len = len(sample)
    # If no negative error, assume errors are symmetric
    if sample_errm is None:
        sample_errm = sample_errp

    # Check if limits
    if sample_ll is None:
        sample_ll = np.zeros(sample_len)
    else:
        if len(sample_ll) != sample_len:
            raise IOError("sample_ll and sample must have same length")
        if ll_max_val is None:
            ll_max_val = (sample.max() + 5 * sample_errm.max()) * np.ones(sample_len)

    if sample_ul is None:
        sample_ul = np.zeros(sample_len)
    else:
        if len(sample_ul) != sample_len:
            raise IOError("sample_ul and sample must have same length")
        if ul_min_val is None:
            ul_min_val = (sample.min() - 5 * sample_errm.max()) * np.ones(sample_len)

    # Create percentiles:
    lower_percentile = (1. - confidence/100.)/2.
    upper_percentile = 1. - lower_percentile

    # Create bins
    if bins is None:
        bin_min = sample.min() - 5 * sample_errm.max()
        bin_max = sample.max() + 5 * sample_errp.max()
        bins = np.linspace(bin_min, bin_max, precision+1)
    else:
        precision = len(bins) - 1
        bin_min = bins.min()
        bin_max = bins.max()

    # Create middle value of the bins, warning this has a length of bins - 1 (i.e. of length precision here)
    # This is essentially used for plotting purposes
    bins_mid = 0.5*(bins[1:]+bins[:-1])

    # Create a array where each line is one realization of MC drawings for each galaxy
    sample_real = np.zeros((N_MC,sample_len))

    # Create array for the weights
    if weights is not None:
        if (len(weights) != sample_len):
            raise IOError("weights and sample must have same length")
        # Create arrays for the realizations of the weights if they have errors
        if weights_err is not None:
            weights_real = np.zeros(sample_real.shape)
            # If list of length 2 provided, interpret as follows asymmetric errors
            if isinstance(weights_err, list) & len(weights_err) == 2:
                weights_errp = weights_err[0]
                weights_errm = weights_err[1]
            # Otherwise symmetric errors
            elif isinstance(weights_err, np.ndarray):
                weights_errp = weights_err
                weights_errm = weights_err
        # If no error on the weights
        else:
            weights_real = weights * np.ones(sample_real.shape)
    # If no weights, use weights = 1 for everything
    else:
        weights_real = np.ones(sample_real.shape)

    # Create array that will hold the PDF for each point in the sample
    sample_pdf = np.zeros((precision_pdf,sample_len))

    if verbose:
        log.info("In compute_CDF_bounds_by_MC: starting Monte Carlo drawings...")
    # For every point in the sample:
    # - Generate N_MC realizations of its value following an asymmetric gaussian pdf with standard deviations being the errors on the point
    # - Calculate its PDF
    sample_real = MC_realization(sample, sample_errp,
                                 sample_errm=sample_errm,
                                 sample_ll=sample_ll,
                                 sample_ul=sample_ul,
                                 ll_max_val=ll_max_val,
                                 ul_min_val=ul_min_val,
                                 N_MC=N_MC,
                                 positive=positive)
    # - If there are weights with errors, generate N_MC realizations of the weight value following same procedure as above
    if weqs:
        weights_real = sample_real.copy()
    elif weights_err is not None and weights is not None:
        for i in range(sample_len):
            weights_real[:,i] = asym_gaussian_draw(weights[i],
                                                   sigma1=weights_errm[i],
                                                   sigma2=weights_errp[i],
                                                   nb_draws=N_MC,
                                                   positive=True)

    # Visualize the PDFs of each individual point in the sample
    if show_plot and (ax is None):
        for i in range(sample_len):
            if (sample_ll[i] == 0) & (sample_ul[i] == 0):
                x, sample_pdf[:,i] = asym_gaussian_pdf(sample[i],
                                                       x_min=bin_min,
                                                       x_max=bin_max,
                                                       precision=len(sample_pdf),
                                                       sigma1=sample_errm[i],
                                                       sigma2=sample_errp[i],
                                                       nb_draws=N_MC,
                                                       positive=positive)
            # If limits, create flat prior over the desired range
            elif (sample_ll[i] == 1) & (sample_ul[i] == 0):
                x = np.linspace(bin_min, bin_max, len(sample_pdf))
                a = 1./np.abs(sample[i] - ll_max_val[i])
                sample_pdf[:,i] = np.zeros(len(sample_pdf))
                sample_pdf[:,i][np.where((x >= sample[i]) & (x <= ll_max_val[i]))] = a
            elif (sample_ll[i] == 0) & (sample_ul[i] == 1):
                x = np.linspace(bin_min, bin_max, len(sample_pdf))
                a = 1./np.abs(ul_min_val[i] - sample[i])
                sample_pdf[:,i] = np.zeros(len(sample_pdf))
                sample_pdf[:,i][np.where((x <= sample[i]) & (x >= ul_min_val[i]))] = a
                if positive:
                    sample_pdf[:,i][np.where(x <= 0)] = 0.0
            else:
                x = np.linspace(bin_min, bin_max, len(sample_pdf))
                a = 1./np.abs(ul_min_val[i] - ll_max_val[i])
                sample_pdf[:,i] = np.zeros(len(sample_pdf))
                sample_pdf[:,i][np.where((x >= ul_min_val[i]) & (x <= ll_max_val[i]))] = a
            ax2.plot(x, sample_pdf[:,i])

        # Create the total PDF and CDF from the individual PDFs of the points in the sample
        summed_pdf = np.sum(sample_pdf, axis=1)
        summed_pdf = summed_pdf/np.sum(summed_pdf*(x[1]-x[0]))
        ax2.plot(x, summed_pdf, c='k', lw=2, label='Summed PDF')
        ax2.legend()

    # Perform bootstrap
    if bootstrap:
        rng = check_random_state(random_state)
        ind = rng.randint(sample_len, size=(N_MC, sample_len))
        sample_bootstrapped = np.zeros(sample_real.shape)
        weights_bootstrapped = np.zeros(weights_real.shape)
        if verbose:
            log.info("In compute_CDF_bounds_by_MC: starting bootstraps...")
        for i in range(N_MC):
            sample_bootstrapped[i] = sample_real[i][ind[i]]
            weights_bootstrapped[i] = weights_real[i][ind[i]]
    else:
        sample_bootstrapped = sample_real
        weights_bootstrapped = weights_real

    if verbose:
        log.info("In compute_CDF_bounds_by_MC: computing CDF...")

    # Compute the PDF and CDF for each realization
    CDF_real = np.zeros((N_MC,precision))
    PDF_real = np.zeros((N_MC,precision))
    for i in range(N_MC):
        hist, bins_ = np.histogram(sample_bootstrapped[i], bins=bins, weights=weights_bootstrapped[i])
        PDF_real[i,:] = hist/float(np.sum(hist)*(bins_[1]-bins_[0]))
        CDF_real[i,:] = np.cumsum(hist).astype(float)/float(np.sum(hist))

    # Compute the percentiles for each bin
    lower = np.zeros(precision)
    median = np.zeros(precision)
    upper = np.zeros(precision)
    for i in range(precision):
        q = mstats.mquantiles(CDF_real[:,i], prob=[lower_percentile, 0.5, upper_percentile])
        median[i] = q[1]
        upper[i] = q[2]
        lower[i] = q[0]

    if show_plot:
        if verbose:
            log.info("In compute_CDF_bounds_by_MC: plotting...")
        if show_median:
            artist, = ax1.plot(bins_mid, median, drawstyle='steps-mid', c=color, **kwargs)
            # Add the first and last line to make the plot look better
            plot_color = plt.getp(artist, 'color')
            plot_ls = plt.getp(artist, 'linestyle')
            plot_lw = plt.getp(artist, 'linewidth')
            plot_zorder = plt.getp(artist, 'zorder')
            imin = np.where(median > 0)[0][0]
            imax = np.where(median == 1)[0][0]
            xmin = bins_mid[imin]
            xmax = bins_mid[imax]
            bottom_line = Line2D([xmin,xmin], [0,median[imin]],
                                 color=plot_color, linestyle=plot_ls, linewidth=plot_lw, zorder=plot_zorder)
            ax1.add_line(bottom_line)
            top_line = Line2D([xmax, xmax+(xmax-xmin)], [median[imax],median[imax]],
                              color=plot_color, linestyle=plot_ls, linewidth=plot_lw, zorder=plot_zorder)
            ax1.add_line(top_line)
        ax1.fill_between(bins_mid, lower, upper, step='mid', color=color, alpha=0.3, zorder=plot_zorder-1)
        ax1.plot(bins_mid, lower, drawstyle='steps-mid',lw=0.7, c=color, zorder=plot_zorder)
        ax1.plot(bins_mid, upper, drawstyle='steps-mid',lw=0.7, c=color, zorder=plot_zorder)

    return bins_mid, median, lower, upper, fig


def MC_realization(sample, sample_errp, sample_errm=None, sample_ll=None, sample_ul=None, ll_max_val=None, ul_min_val=None, N_MC=1000, positive=False):
    """
        Function to create a realization of a sample with errors and upper limits.
        This assumes the value's PDF's can be represented by asymmetric gaussians whose sigmas are the plus and minus error.
        For limits, assume a flat prior (uniform draw) between the limit and the array specified by ll_max_val/ul_min_val.
    """

    # Check if limits
    if sample_ll is None:
        sample_ll = np.zeros(len(sample))
    else:
        if len(sample_ll) != len(sample):
            raise IOError("sample_ll and sample must have same length")
        if ll_max_val is None:
            ll_max_val = (sample.max() + 5 * sample_errm.max()) * np.ones(len(sample))

    if sample_ul is None:
        sample_ul = np.zeros(len(sample))
    else:
        if len(sample_ul) != len(sample):
            raise IOError("sample_ul and sample must have same length")
        if ul_min_val is None:
            ul_min_val = (sample.min() - 5 * sample_errm.max()) * np.ones(len(sample))

    sample_real = np.zeros((N_MC,len(sample)))
    for i in range(len(sample)):
        # If no limits draw asymmetric gaussian
        if (sample_ll[i] == 0) & (sample_ul[i] == 0):
            sample_real[:,i] = asym_gaussian_draw(sample[i],
                                                  sigma1=sample_errm[i],
                                                  sigma2=sample_errp[i],
                                                  nb_draws=N_MC,
                                                  positive=positive)
        # Otherwise, draw uniform
        elif (sample_ll[i] == 1) & (sample_ul[i] == 0):
            sample_real[:,i] = np.random.rand(N_MC) * (ll_max_val[i] - sample[i]) + sample[i]
        elif (sample_ll[i] == 0) & (sample_ul[i] == 1):
            if positive & (ul_min_val[i] < 0.0):
                sample_real[:,i] = np.random.rand(N_MC) * sample[i]
            else:
                sample_real[:,i] = np.random.rand(N_MC) * (sample[i] - ul_min_val[i]) + ul_min_val[i]
        else:
            sample_real[:,i] = np.random.rand(N_MC) * (ll_max_val[i] - ul_min_val[i]) + ul_min_val[i]

    return sample_real


def asym_gaussian_draw(mu, sigma1, sigma2, nb_draws=1000, precision=500, positive=False, ax=None, **kwargs):
    """
    Function that draws randomly in a asymmetric gaussian distribution.
    in the form :   { exp( -(x-mu)**2/(2*sigma1**2) )     if x < mu
                    { exp( -(x-mu)**2/(2*sigma2**2) )     if x >= mu
    Also plots the distribution if ax is not None
    Returns an array of the drawings
    """

    if sigma1 == 0. and sigma2 == 0.:
        return mu * np.ones(nb_draws)

    if sigma1 < 0. or sigma2 < 0.:
        raise ValueError('sigma1 or sigma2 can not be negative, check your input.')

    if sigma1 == 0.:
        if mu == 0:
            sigma1 = 1e-9
            log.warning('In asym_gaussian_draw: sigma1 and mu are equal to zero, replacing sigma1 by 1e9')
        else:
            sigma1 = 1e-9 * mu
            log.warning('In asym_gaussian_draw: sigma1 is equal to zero, replacing sigma1 by mu * 1e9')
    if sigma2 == 0.:
        if mu == 0:
            sigma2 = 1e-9
            log.warning('In asym_gaussian_draw: sigma2 and mu are equal to zero, replacing sigma2 by 1e9')
        else:
            sigma2 = 1e-9 * mu
            log.warning('In asym_gaussian_draw: sigma2 is equal to zero, replacing sigma2 by mu * 1e9')

    # limits
    x_min = mu - 10.*sigma1
    x_max = mu + 10.*sigma2

    # create Cumulative distribution
    x = np.linspace(x_min, x_max, precision)
    p_x = anpdf(x, mu, sigma1, sigma2)
    F_x = ancdf(x, mu, sigma1, sigma2)

    # draw from generated distribution
    if positive:
        min_rand = F_x[x.searchsorted(0)]
        rand = (1-min_rand) * np.random.rand(nb_draws) + min_rand
    else:
        rand = np.random.rand(nb_draws)
    draw = x[F_x.searchsorted(rand)]

    if ax is not None:
        label = '\n'.join([r'$\mu = ${:.2f}'.format(mu),
                           r'$\sigma_- = ${:.2f}'.format(sigma1),
                           r'$\sigma_+ = ${:.2f}'.format(sigma2)])
        ax.hist(draw, bins=20, density=True, label=label, **kwargs)
        ax.plot(x, p_x, color='k')
        ax.legend()

    return draw


def asym_gaussian_pdf(mu, sigma1, sigma2, x_min=None, x_max=None, cumulative=False, precision=500, ax=None, **kwargs):
    """
    Function that returns an asymmetric gaussian distribution.
    in the form :   { exp( -(x-mu)**2/(2*sigma1**2) )     if x < mu
                    { exp( -(x-mu)**2/(2*sigma2**2) )     if x >= mu
    Also plots the distribution if ax is not None
    Returns an array of the probability distribution function and its abscissa.
    """
    # Normalization
    # norm = 1. / ( np.sqrt(np.pi/2.) * (sigma1 + sigma2) )

    if sigma1 < 0. or sigma2 < 0.:
        raise ValueError('sigma1 or sigma2 can not be negative, check your input.')

    if sigma1 == 0.:
        if mu == 0:
            sigma1 = 1e-9
            log.warning('In asym_gaussian_pdf: sigma1 and mu are equal to zero, replacing sigma1 by 1e9')
        else:
            sigma1 = 1e-9 * mu
            log.warning('In asym_gaussian_pdf: sigma1 is equal to zero, replacing sigma1 by mu * 1e9')
    if sigma2 == 0.:
        if mu == 0:
            sigma2 = 1e-9
            log.warning('In asym_gaussian_pdf: sigma2 and mu are equal to zero, replacing sigma2 by 1e9')
        else:
            sigma2 = 1e-9 * mu
            log.warning('In asym_gaussian_pdf: sigma2 is equal to zero, replacing sigma2 by mu * 1e9')

    # limits
    if x_min is None:
        x_min = mu - 10.*sigma1
    if x_max is None:
        x_max = mu + 10.*sigma2

    # create Cumulative distribution
    x = np.linspace(x_min, x_max, precision)
    F_x = ancdf(x, mu=mu, s1=sigma1, s2=sigma2)
    p_x = anpdf(x, mu=mu, s1=sigma1, s2=sigma2)

    if ax is not None:
        ax.plot(x, p_x, **kwargs)
        ax.legend()

    output = p_x
    if cumulative:
        output = F_x

    return x, output


def gaussian(x, mu, sigma, normed=False):
    if sigma <= 0:
        raise ValueError("invalid value for sigma in gaussian, must be non-zero.")
    gaussian = np.exp(-(x-mu)**2/(2*sigma**2))
    if normed:
        norm = 1./np.sqrt(2*np.pi*sigma**2)
    else:
        norm = 1.
    return norm*gaussian


def anpdf(x, mu=0, s1=1, s2=1):
    a = 2. / (np.sqrt(2. * np.pi) * (s1 + s2))
    if isinstance(x, np.ndarray):
        pdf = a * np.exp(- 0.5 * (x - mu) ** 2 / s1 ** 2)
        w = np.where(x > mu)
        pdf[w] = a * np.exp(- 0.5 * (x[w] - mu) ** 2 / s2 ** 2)
    else:
        if x < mu:
            pdf = a * np.exp(- 0.5 * (x - mu) ** 2 / s1 ** 2)
        else:
            pdf = a * np.exp(- 0.5 * (x - mu) ** 2 / s2 ** 2)
    return pdf


def ancdf(x, mu=0, s1=1, s2=1):
    from scipy.stats import norm
    a = 2. / (s1 + s2)
    if isinstance(x, np.ndarray):
        cdf = a * s1 * norm.cdf(x, mu, s1)
        w = np.where(x > mu)
        cdf[w] = a * (0.5 * (s1 - s2) + s2 * norm.cdf(x[w], mu, s2))
    else:
        if x < mu:
            cdf = a * s1 * norm.cdf(x, mu, s1)
        else:
            cdf = a * (0.5 * (s1 - s2) + s2 * norm.cdf(x, mu, s2))
    return cdf
