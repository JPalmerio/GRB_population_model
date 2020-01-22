import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotting_functions as pf
from io_grb_pop import read_column, root_dir
from miscelaneous import chi2_func, log_to_lin
from stats import unbinned_empirical_cdf
import logging

log = logging.getLogger(__name__)
T_live_BATSE = 6.54  # years (Antier-Farfar thesis ?)
T_mission = 9.1  # years (Stern et al. 2001)
Omega_BATSE = 0.67 * 4. * np.pi  # years (Antier-Farfar thesis ?)


def double_array(array_in):
    """
        Convenience function to duplicate the consecutive rows in an
        array to make it easy to plot.
    """
    array_out = np.zeros(2*len(array_in))
    for i in range(len(array_in)):
        array_out[2*i] = array_in[i]
        array_out[2*i+1] = array_in[i]
    return array_out


def plottable_hist(x_in, y_in, last_bin_edge):
    """
        Convenience function to create a plottable pair of x,y arrays
        - x is the left edge of the bins so last_bin_edge must be
        provided
        - y[i] is the histogram value between x[i] and x[i+1]
        The output can be easily plotted.
        Example:
        x = [0,1,2]; last_bin_edge = 3
        y = [2,4,3]
        ----> x_out = [0,1,1,2,2,3] y_out = [2,2,4,4,3,3]
    """
    if len(x_in) != len(y_in):
        print("[ERROR] in plottable_hist: x and y have different lenths : {} and {}".format(len(x_in), len(y_in)))
        raise IOError

    x_out = double_array(x_in)
    y_out = double_array(y_in)

    x_out[1:-1] = x_out[2:]
    x_out[0] = x_in[0]
    x_out[-1] = last_bin_edge

    return x_out, y_out


def create_Stern_hist(fname=None, verbose=False):
    """
        Create the histogram from [Stern et al. 2001]
        (https://ui.adsabs.harvard.edu/abs/2001ApJ...563...80S/abstract)
        Fig. 23 rebinned to avoid low counts
    """
    global T_live_BATSE
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


def read_SHOALS_file(SHOALS_file=None):
    if SHOALS_file is None:
        SHOALS_file = root_dir/'catalogs/SHOALS_cat/SHOALS_cat.txt'

    df_obs = pd.read_csv(SHOALS_file, sep='\t', header=3, low_memory=False)
    keys = ['S_BAT', 'z']
    for key in keys:
        df_obs[key] = pd.to_numeric(df_obs[key], errors='coerce')
    df_obs = df_obs.dropna()
    return df_obs


def plot_constraint_old(constraint, bins, mod, obs, err):
    """
        Convenience function to quickly plot the various constraints
        used in the population model
        The histograms are doubled for plotting (easier...)
    """
    if constraint not in ['Stern', 'EpGBM', 'eBAT6']:
        raise ValueError("constraint must be one of {}".format(['Stern', 'EpGBM', 'eBAT6']))

    x, mod_to_plot = plottable_hist(bins[:-1], mod, last_bin_edge=bins[-1])
    x, obs_to_plot = plottable_hist(bins[:-1], obs, last_bin_edge=bins[-1])
    err_to_plot = double_array(err)

    fig, ax = plt.subplots()

    ax.plot(x, mod_to_plot, drawstyle='steps-post', label='Model')
    ax.plot(x, obs_to_plot, drawstyle='steps-post', label=constraint, color='k')
    ax.fill_between(x, obs_to_plot-err_to_plot, obs_to_plot+err_to_plot,
                    step='post', color='k', alpha=0.3)
    if constraint == 'Stern':
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel(r'Peak flux $\rm{[ph\,cm^{-2}\,s^{-1}\,50-300\,keV]}$')
        ax.set_ylabel(r'$\Delta N/\Delta \,(\rm{log\,}P)\,\rm{[yr^{-1}] }$ in $4\pi$')
    elif constraint == 'EpGBM':
        ax.set_xlabel(r'Peak Energy $\rm{[keV]}$')
        ax.set_ylabel(r'Number of GRBs')
        ax.set_xlim([13.,10000.])
        ax.set_ylim([0.,350.])
        ax.set_xscale('log')
    ax.legend()
    return


def plot_EpGBM(fname=None, log=True, density=False, ax=None, **kwargs):
    if fname is None:
        fname = root_dir/'observational_constraints/Ep_GBM_for_plotting.txt'
    if ax is None:
        ax = plt.gca()
    bins, obs_lin, err_lin = create_EpGBM_hist(fname, density=density)

    x, x_errp, x_errm = xerr_from_bins(bins, log=log)

    ax.errorbar(x, obs_lin, xerr=[x_errm, x_errp], yerr=err_lin, **kwargs)

    return


def xerr_from_bins(bins, log=False):
    """
        Convenience function to easily calculate the x_err from bins to
        represent a histogram as crosses on a plot
        returns:
        - x the middle point of the bin in the x direction
        - x_errp the extent of the bin from x up to the right bin edge
        - x_errm the extent of the bin from x down to the left bin edge
    """
    if log:
        x = np.sqrt(bins[1:] * bins[:-1])
    else:
        x = 0.5 * (bins[1:] + bins[:-1])

    x_errp = bins[1:] - x
    x_errm = x - bins[:-1]

    return x, x_errp, x_errm


def plot_constraint(constraint, bins, mod, obs, err, plot_mod=True, plot_obs=True, bins_log=True, ax=None):
    """
        Convenience function to quickly plot the various constraints
        used in the population model
        The histograms are doubled for plotting (easier...)
    """
    if constraint not in ['Stern', 'EpGBM', 'eBAT6']:
        raise ValueError("constraint must be one of {}".format(['Stern', 'EpGBM', 'eBAT6']))

    x, x_errp, x_errm = xerr_from_bins(bins, log=log)

    if ax is None:
        fig, ax = plt.subplots()

    if plot_mod:
        x1, mod_to_plot = plottable_hist(bins[:-1], mod, last_bin_edge=bins[-1])
        ax.plot(x1, mod_to_plot, drawstyle='steps-post', label='Model', lw=2)

    if plot_obs:
        x1, obs_to_plot = plottable_hist(bins[:-1], obs, last_bin_edge=bins[-1])
        ax.errorbar(x, obs, xerr=[x_errm, x_errp], yerr=err,
                    capsize=0, label='Observations', fmt='.', color='k')

    if constraint == 'Stern':
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel(r'Peak flux $\rm{[ph\,cm^{-2}\,s^{-1}\,50-300\,keV]}$')
        ax.set_ylabel(r'$\Delta N/\Delta \,(\rm{log\,}P)\,\rm{[yr^{-1}] }$ in $4\pi$')
    elif constraint == 'EpGBM':
        ax.set_xlabel(r'Peak Energy $\rm{[keV]}$')
        ax.set_ylabel(r'Number of GRBs')
        ax.set_xlim([13.,10000.])
        ax.set_ylim([0.,350.])
        ax.set_xscale('log')
    elif constraint == 'eBAT6':
        ax.set_xlabel('Redshift(z)')
        ax.set_ylabel(r'Number of GRBs')
        ax.set_xlim([0.01, 6.])
    ax.legend()
    return


def create_EpGBM_hist(fname=None, verbose=False, density=False):
    """
        Create the histogram from Gruber et al. 2014 catalog of GBM
        bursts with a pflx cut of 0.9 ph/s/cm2
        If density is true, hist is returned as :
        delta_N / (N_EpGBM * delta_log_bin)
    """
    if fname is None:
        fname = root_dir/'observational_constraints/Ep_GBM_for_plotting.txt'

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


def create_eBAT6_hist(fname, verbose=False, eBAT6_weight=10):
    """
        Create the histogram for the eBAT6 redshift distribution
        Data was taken from Pescalli et al. 2016
    """

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


def plot_SHOALS_flnc_z_plane(df_mod, SHOALS_file=None, bins_flnc=None, bins_z=None):
    """
        Compare the model predictions with the observed SHOALS sample in the erg fluence vs redshift
        plane. Sample is defined as GRBs with energy fluence in 15-150 keV above 1e-6 erg/cm2
        Data is compiled from Perley+16
    """

    df_obs = read_SHOALS_file(SHOALS_file)

    if bins_flnc is None:
        bins_flnc = np.linspace(-6, -3.5, 15)
    if bins_z is None:
        bins_z = np.linspace(0,7, 14)

    fig, axes = pf.cool_hist2d(x=df_mod['z'].to_numpy(),
                               y=np.log10(df_mod['erg_flnc_BAT'].to_numpy()), cb=False,
                               mode='hist2d',
                               left_hist_kwargs={'color': 'lightgray','alpha':0.6, 'bins':bins_flnc},
                               top_hist_kwargs={'color': 'lightgray','alpha':0.6, 'bins':bins_z},
                               top_kdeplot_kwargs={'color': 'k', 'label': 'Model'},
                               left_kdeplot_kwargs={'color': 'k', 'label': None})

    pf.cool_hist2d(x=df_obs['z'].to_numpy(),
                   y=np.log10(1e-7*df_obs['S_BAT'].to_numpy()),
                   fig=fig, cb=False,
                   mode='scatter',
                   left_kdeplot_kwargs={'color': 'C5', 'label': None},
                   left_hist_kwargs={'color': 'C5','alpha':0.5, 'bins':bins_flnc},
                   top_hist_kwargs={'color': 'C5','alpha':0.5, 'bins':bins_z},
                   top_kdeplot_kwargs={'color': 'C5', 'label': 'SHOALS observed'},
                   color='C5', linewidth=0.5, label='SHOALS observed')

    axes['center'].set_ylim(bins_flnc[0], bins_flnc[-1])
    axes['center'].set_xlim(bins_z[0], bins_z[-1])
    axes['center'].set_xlabel('Redshift (z)')
    axes['left'].set_ylabel('log erg fluence [erg/cm2]')
    axes['center'].legend()

    return


def plot_observed_SHOALS_z_distr(ax, SHOALS_file=None, **kwargs):
    df_obs = read_SHOALS_file(SHOALS_file)
    z_obs, cdf_obs = unbinned_empirical_cdf(df_obs['z'].to_numpy())
    ax.plot(z_obs, cdf_obs, drawstyle='steps-post', **kwargs)
    return


def plot_SHOALS_distr(df_mod, key, SHOALS_file=None, cumul=False, fig=None, plot_obs=True,
                      mod_color='C0', mod_label='Model', bins=None):

    if fig is None:
        fig, ax = plt.subplots(figsize=(10,8))
    else:
        ax = fig.axes[0]

    df_obs = read_SHOALS_file(SHOALS_file)
    if cumul:
        if plot_obs:
            if key == 'S_BAT':
                data_obs = np.log10(df_obs[key].to_numpy())
            else:
                data_obs = df_obs[key].to_numpy()

            x_obs, cdf_obs = unbinned_empirical_cdf(data_obs)
            ax.plot(x_obs, cdf_obs, color='C5', drawstyle='steps-post', label='SHOALS observed', lw=3)

        x_mod, cdf_mod = unbinned_empirical_cdf(df_mod.to_numpy())
        ax.plot(x_mod, cdf_mod, color=mod_color, drawstyle='steps-post', label=mod_label)

        ax.legend()
        ax.set_ylabel('CDF')
        ax.set_ylim(0,1)
    else:
        if bins is None:
            bins = np.linspace(0,7,14)
        if plot_obs:
            sns.kdeplot(df_obs[key], ax=ax, label='SHOALS observed', color='C5')
            ax.hist(df_obs[key], color='C5', bins=bins, density=True, alpha=0.7)

        sns.kdeplot(df_mod, ax=ax, label='Model', color=mod_color)
        ax.hist(df_mod, color=mod_color, bins=bins, density=True, alpha=0.7)

        ax.set_ylabel('PDF')
    if key == 'z':
        ax.set_xlabel('Redshift (z)')
        ax.set_xlim(0,8)
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


def compare_to_Stern(Stern_sample, Stern_file, method='chi2',
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
        plot_constraint('Stern', bins, norm*mod, obs_lin, err_lin)
    return norm, norm*mod


def compare_to_EpGBM(GBM_sample, GBM_file, method='chi2', show_plot=False):
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
        plot_constraint('EpGBM', bins, norm*mod, obs_lin, err_lin)

    return chi2


def compare_to_eBAT6(eBAT6_sample, eBAT6_file, method='chi2', show_plot=False):
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
        plot_constraint('eBAT6', bins, norm*mod, obs_lin, err_lin, bins_log=False)

    return chi2
