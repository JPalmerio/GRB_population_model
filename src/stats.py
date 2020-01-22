import time
import logging
import numpy as np
import pandas as pd
import physics as ph
import plotting_functions as pf
from scipy.stats import mstats, ks_2samp
from GRB_population import GRBPopulation
# from observational_constraints import compare_to_Stern

log = logging.getLogger(__name__)


def MonteCarlo_routine(Nb_GRBs, cosmo, params, incl_samples, incl_instruments, obs_dir,
                       ECLAIRs_prop=None, output_dir=None, verbose=False, run_mode=None):

    log.debug("Starting drawings of random samples...")
    t1 = time.time()
    GRB_population = GRBPopulation(Nb_GRBs, output_dir=output_dir)
    GRB_prop = GRB_population.draw_GRB_properties(cosmo=cosmo, params=params, run_mode=run_mode)
    t2 = time.time()
    log.debug(f"Drawings done in {t2-t1:.3f} s")

    ph.calc_peak_photon_flux(GRB_prop, incl_instruments, ECLAIRs_prop)
    ph.calc_peak_energy_flux(GRB_prop, incl_instruments, ECLAIRs_prop)
    ph.calc_photon_fluence(GRB_prop, incl_instruments)
    ph.calc_energy_fluence(GRB_prop, incl_instruments)
    ph.calc_det_prob(GRB_prop, incl_samples, **ECLAIRs_prop)

    df = pd.DataFrame(GRB_prop)

    # if 'Stern' in incl_samples:
    #     norm = compare_to_Stern(df[(df['pdet_Stern'] == 1)]['pht_pflx_BATSE'],
    #                             Stern_file=obs_dir/'Stern_lognlogp_rebinned.txt',
    #                             Nb_GRBs=Nb_GRBs, show_plot=verbose)
    return


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


def subsample_and_KS(df1, df2, N_sub, key, confidence=95.0, N_bs=100, precision=500, bins=None, show_plot=False,
                     label1=None, label2=None, subsample1=True, subsample2=True):
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

    t1 = time.time()
    CDF1 = CDF_with_bootstrap(subsamp1, bins=bins)
    CDF2 = CDF_with_bootstrap(subsamp2, bins=bins)
    med1, lw1, up1 = compute_CDF_quantiles(CDF1, confidence=confidence)
    med2, lw2, up2 = compute_CDF_quantiles(CDF2, confidence=confidence)
    t2 = time.time()
    log.info(f"CDF calculation done in {t2-t1:.3f} s")

    pfrac = len(p_value[np.where(p_value < (1.-confidence/100.))[0]])/len(p_value)

    if show_plot:
        pf.plot_CDFs_and_KS_results(bins, med1, med2, lw1, lw2, up1, up2, D_stat, p_value, confidence,
                                    pfrac, label1, label2)

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


def CDF_with_bootstrap(subsample, bins):
    """
        Takes in a subsample array of size (N_bs, N_sub) and computes the CDF for each bootstrap realization
        The CDFs are sampled given a certain precision determined from the bins.
        N_bs stands for number of bootstraps
        N_sub stands for size of the subsample
    """
    N_bs = subsample.shape[0]
    precision = bins.shape[0] - 1
    # Compute the CDF for each subsample realization
    CDF_real = np.zeros((N_bs, precision))
    for i in range(N_bs):
        hist, bins_ = np.histogram(subsample[i], bins=bins)
        CDF_real[i,:] = np.cumsum(hist).astype(float)/float(np.sum(hist))

    return CDF_real


def compute_CDF_quantiles(CDF, confidence=95.0):
    """
        Takes a 2D array of CDFs of size (N_bs, N_bins).
        N_bs stands for the number of bootstraps
        N_bins stands for the number of bins within the CDF.
        Returns the median, lower and upper bounds at the desired confidence level.
    """

    N_bins = CDF.shape[1]
    # Create percentiles:
    lower_percentile = (1. - confidence/100.)/2.
    upper_percentile = 1. - lower_percentile
    # Compute the percentiles for each bin
    lower = np.zeros(N_bins)
    median = np.zeros(N_bins)
    upper = np.zeros(N_bins)
    for i in range(N_bins):
        q = mstats.mquantiles(CDF[:,i], prob=[lower_percentile, 0.5, upper_percentile])
        median[i] = q[1]
        upper[i] = q[2]
        lower[i] = q[0]
    return median, lower, upper


def compute_KS_with_bootstrap(sample1, sample2, N_bs):
    """
        Compute the 2 sample K-S test N_bs times where
        N_bs stands for the number of bootstrap samples
        Return the D-statistic and the pvalue for each bootstrapped
        realization in the form of 1D arrays of size N_bs
    """
    D_stat = np.zeros(N_bs)
    p_value = np.zeros(N_bs)
    for i in range(N_bs):
        D_stat[i], p_value[i] = ks_2samp(sample1[i], sample2[i])
    return D_stat, p_value
