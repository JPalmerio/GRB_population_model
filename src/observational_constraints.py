import numpy as np
import plotting_functions as pf
from io_grb_pop import read_column, root_dir
from miscelaneous import chi2_func, log_to_lin
import logging

log = logging.getLogger(__name__)
T_live_BATSE = 6.54  # years (Antier-Farfar thesis ?)
T_mission = 9.1  # years (Stern et al. 2001)
Omega_BATSE = 0.67 * 4. * np.pi  # years (Antier-Farfar thesis ?)


def create_Stern_hist(fname=None, verbose=False):
    """
        Create the histogram from [Stern et al. 2001]
        (https://ui.adsabs.harvard.edu/abs/2001ApJ...563...80S/abstract)
        Fig. 23 rebinned to avoid low counts
    """
    global T_live_BATSE
    if fname is None:
        fname = root_dir/'observational_constraints/Stern_lognlogp_rebinned.txt'

    bins = read_column(fname, 0, array=False)
    hist = read_column(fname, 1)
    err = read_column(fname, 2)
    bins.append(50.0)  # append right edge of last bin
    bins = np.array(bins)

    if verbose:
        ln_oi = 0.
        for i, val in enumerate(hist):
            val2 = 10**val
            delta_bin = np.log10(bins[i+1]/bins[i])
            val2 *= delta_bin * T_live_BATSE
            ln_oi += val2*np.log(val2) - val2
        print(f"ln(o_i!) = {ln_oi} from Stern histogram")

    return bins, hist, err


def create_EpGBM_hist(fname=None, verbose=False, density=False):
    """
        Create the histogram from Gruber et al. 2014 catalog of GBM
        bursts with a pflx cut of 0.9 ph/s/cm2
        If density is true, hist is returned as :
        delta_N / (N_EpGBM * delta_log_bin)
    """
    if fname is None:
        fname = root_dir/'observational_constraints/EpGBM_for_plotting.txt'

    bins = read_column(fname, 0, array=False)
    hist = read_column(fname, 1)
    err = read_column(fname, 2)
    bins.append(1e4)  # append right edge of last bin
    bins = np.array(bins)

    if density:
        N_EpGBM = hist.sum()
        delta_bin = np.log10(bins[1:]/bins[:-1])
        hist /= (N_EpGBM * delta_bin)
    if verbose:
        ln_oi = 0.
        for i, val in enumerate(hist):
            ln_oi += val*np.log(val) - val
        print(f"ln(o_i!) = {ln_oi} from EpGBM histogram")

    return bins, hist, err


def create_eBAT6_hist(fname=None, verbose=False, eBAT6_weight=10):
    """
        Create the histogram for the eBAT6 redshift distribution
        Data was taken from Pescalli et al. 2016
    """
    if fname is None:
        fname = root_dir/'observational_constraints/eBAT6_constraint_for_plotting.txt'

    bins = read_column(fname, 0, array=False)
    hist = read_column(fname, 1)
    err = read_column(fname, 2)
    bins.append(6)  # append right edge of last bin
    bins = np.array(bins)
    if verbose:
        ln_oi = 0.
        for i, val in enumerate(hist):
            if val != 0:
                ln_oi += val*np.log(val) - val
        ln_oi *= eBAT6_weight
        print(f"ln(o_i!) = {ln_oi} from eBAT6 histogram")

    return bins, hist, err


def plot_EpL_plane(GRB_pop, sample='eBAT6'):
    """
        Plot the GRB population in the Ep-L plane
    """
    fig, axes = pf.fig_marg()
    mini_cbax = fig.add_axes([0.780, 0.125, 0.015, 0.20])
    raise NotImplementedError
    return


def global_GRB_rate_Stern(Stern_file):
    """
        Calculate the global LGRB rate by correcting the Stern histogram
        with the efficiency correction [Stern et al. 2001]
        (https://ui.adsabs.harvard.edu/abs/2001ApJ...563...80S/abstract)
        Eq. 5. and summing over all the bins
    """
    global T_live_BATSE, Omega_BATSE
    bins, hist_obs, err_obs = create_Stern_hist(Stern_file)
    bin_midpoint = 0.5*(bins[:-1]+bins[1:])
    hist_obs_corrected = hist_obs / efficiency_correction_Stern(bin_midpoint)

    N_GRB_BATSE_tot = np.sum(hist_obs_corrected)
    glob_rate = N_GRB_BATSE_tot / T_live_BATSE
    all_sky_glob_rate = glob_rate * Omega_BATSE

    log.info(''.join([f"Global LGRB rate from Stern constraint: {all_sky_glob_rate:.2f} ",
             f"GRB/yr in 4 pi with peak flux in [50-300 keV] above {bins.min()} ph/s/cm2"]))

    return all_sky_glob_rate


def efficiency_correction_Stern(pflx, c_e0=0.097, nu=2.34, norm=0.7):
    """
        The efficiency function of BATSE for detecting GRBs as a function of peak flux, derived by Stern+01
        c_e0 is in [counts/s/cm2]
        pflx is in [ph/s/cm2]
    """
    c_e = pflx * 0.75  # the conversion factor from counts to pflx comes from the Stern+01 paper as well, figure 7.
    return norm * (1.0 - np.exp(-(c_e/c_e0)**2))**nu


def compare_to_Stern(Stern_sample, Stern_file=None, method='chi2',
                     show_plot=False, Nb_GRBs=None, new=False, Stern_det_prob=None):
    """
        Compare the predictions of the model with the observations of
        [Stern et al. 2001]
        (https://ui.adsabs.harvard.edu/abs/2001ApJ...563...80S/abstract)
        Fig. 23.
        Calculate the chi2 between the model and the observations.
        Use chi2 minimization to get the normalization of the model.
        Returns the normalization and the value of the Stern histogram
        predicted by the model with this normalization applied to it.
    """
    global T_live_BATSE, Omega_BATSE

    if method != 'chi2':
        raise NotImplementedError("Only Chi2 is implemented for now")

    if new:
        if Stern_det_prob is None:
            raise IOError('If using the new Stern mode, Stern_det_prob must be provided')
        # Assumes the file for Stern has linear data
        bins, obs_lin, err_lin = create_Stern_hist(Stern_file)
        mod, _u = np.histogram(Stern_sample, weights=Stern_det_prob, bins=bins)
        norm = np.sum(obs_lin*mod/err_lin**2) / np.sum((mod/err_lin)**2)
        chi2 = chi2_func(norm*mod, obs_lin, err_lin)
        norm *= Omega_BATSE / T_live_BATSE
    else:
        bins, obs_log, err_log = create_Stern_hist(Stern_file)
        obs_lin, err_lin, _u = log_to_lin(obs_log, err_log)
        delta_bin = np.log10(bins[1:]/bins[:-1])
        grate = np.sum(delta_bin * obs_lin)
        log.info("Global GRB rate: {}".format(grate))
        mod_raw, _u = np.histogram(Stern_sample, bins=bins)
        mod = mod_raw / delta_bin
        norm = np.sum(obs_lin*mod/err_lin**2) / np.sum((mod/err_lin)**2)
        chi2 = chi2_func(norm*mod, obs_lin, err_lin)

    log.info(f"          Stern chi2: {chi2:.4e}")
    log.info(f" Stern normalization: {norm:.4e} yr-1")
    log.info(f" Simulation duration: {1./norm:.4e} yr")
    log.info(f"           ndot_LGRB: {norm*Nb_GRBs:.4e} LGRB/yr")
    if show_plot:
        pf.plot_constraint('Stern', bins, norm*mod, obs_lin, err_lin)
    return norm, norm*mod


def compare_to_EpGBM(GBM_sample, GBM_file=None, method='chi2', show_plot=False):
    """
        Compare the predictions of the model with the EpGBM sample
        defined as GRBs with peak flux in 50-300 keV above 0.9 ph/s/cm2
        Data is compiled from Gruber+14
    """

    if method != 'chi2':
        raise NotImplementedError("Only Chi2 is implemented for now")

    bins, obs_lin, err_lin = create_EpGBM_hist(GBM_file)
    mod, _u = np.histogram(np.log10(GBM_sample), bins=np.log10(bins))
    norm = np.sum(obs_lin*mod/err_lin**2) / np.sum((mod/err_lin)**2)
    chi2 = chi2_func(norm*mod, obs_lin, err_lin)
    log.info(f"          EpGBM chi2: {chi2:.4e}")
    log.info(f" EpGBM normalization: {norm:.4e}")
    if show_plot:
        pf.plot_constraint('EpGBM', bins, norm*mod, obs_lin, err_lin)

    return chi2


def compare_to_eBAT6(eBAT6_sample, eBAT6_file=None, method='chi2', show_plot=False):
    """
        Compare the predictions of the model with the eBAT6 sample
        defined as GRBs with peak flux in 15-150 keV above 2.6 ph/s/cm2
        Data is compiled from Pescalli+16
    """

    if method != 'chi2':
        raise NotImplementedError("Only Chi2 is implemented for now")

    bins, obs_lin, err_lin = create_eBAT6_hist(eBAT6_file)
    mod, _u = np.histogram(eBAT6_sample, bins=bins)
    norm = np.sum(obs_lin*mod/err_lin**2) / np.sum((mod/err_lin)**2)
    chi2 = chi2_func(norm*mod, obs_lin, err_lin)
    log.info(f"          eBAT6 chi2: {chi2:.4e}")
    log.info(f" eBAT6 normalization: {norm:.4e}")
    if show_plot:
        pf.plot_constraint('eBAT6', bins, norm*mod, obs_lin, err_lin, bins_log=False)

    return chi2
