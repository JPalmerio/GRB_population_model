import corner
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import miscellaneous as msc
import stats as st
import io_grb_pop as io
from io_grb_pop import root_dir
from constants import T_live_BATSE


def plot_intensity_constraint(ax, pop=None, plot_obs=True, alt_norm=False, **kwargs):
    """
        Plot the intensity constraint from Stern et al. 2001
    """

    if pop is not None:
        # Turn normal hist into logN-logP
        bins = pop.mock_constraints['Stern']['bins']
        delta_bin = np.log10(bins[1:]/bins[:-1])
        if alt_norm:
            err = pop.mock_constraints['Stern']['err']
            hist = pop.mock_constraints['Stern']['hist']
            hist = hist / (T_live_BATSE*delta_bin)
            err = err / (T_live_BATSE*delta_bin)
        else:
            err = pop.mock_constraints['Stern']['err_unnormed']
            hist = pop.mock_constraints['Stern']['hist_unnormed']
            hist = hist / (pop.normalization['T_sim']*delta_bin)
            err = err / (pop.normalization['T_sim']*delta_bin)

        y_mod, y_mod_errp, y_mod_errm = msc.lin_to_log(hist, err)

        # Make it look better for plots
        x_plottable, y_mod_plottable = plottable_hist(bins[:-1], y_mod, last_bin_edge=bins[-1])
        y_mod_errm_plottable = double_array(y_mod_errm)
        y_mod_errp_plottable = double_array(y_mod_errp)

        # Plot them
        art, = ax.plot(x_plottable, y_mod_plottable, **kwargs)
        ax.fill_between(x_plottable,
                        y_mod_plottable+y_mod_errp_plottable,
                        y_mod_plottable-y_mod_errm_plottable,
                        alpha=0.6, color=plt.getp(art,'color'), zorder=plt.getp(art,'zorder'))
    if plot_obs:
        bins_Stern, obs_log, err_log = io.read_logRlogN()
        x, x_errp, x_errm = xerr_from_bins(bins_Stern, logscale=True)
        ax.errorbar(x, obs_log, xerr=[x_errm, x_errp], yerr=err_log, color='k', fmt='none')

    ax.set_xscale('log')
    ax.set_xlabel(r'$N^{\rm pk}_{50-300\,\rm keV}~\rm{[ph\,cm^{-2}\,s^{-1}]}$')
    ax.set_ylabel(r'log $\Delta R/\Delta \,(\log\,N^{\rm pk})$')
    return


def plot_spectral_constraint(ax, pop=None, plot_obs=True, **kwargs):

    if plot_obs:
        bins_Ep = np.linspace(0,4.5, 30)
        bins, obs_lin, err_lin = io.read_constraint(name='EpGBM', density=True, bins_log=True, last_bin_edge=1e4)
        x, x_errp, x_errm = xerr_from_bins(bins, logscale=True)
        ax.errorbar(x, obs_lin, xerr=[x_errm, x_errp], yerr=err_lin, color='k', fmt='none')
        plot_obs_property(root_dir/'catalogs/GBM_cat/fermi_GBM_cat_total.dat',
                          key='pflx_band_epeak',
                          func=[msc.filter_df, msc.filter_df],
                          func_args=[{'filtering_key':'t90', 'lim_min':2, 'errors':'coerce'},
                                     {'filtering_key':'pflx_band_phtfluxb', 'lim_min':0.9, 'errors':'coerce'}],
                          log=True, verbose=True, kde=True,
                          bins=bins_Ep,
                          ax=ax, density=True, label='Spectral Constraint', alpha=0.8, zorder=0,
                          color='lightgray', errors='coerce')
        ax.get_legend().remove()
    if pop is not None:
        # Read the mock constraint from the population
        bins = pop.mock_constraints['EpGBM']['bins']
        err = pop.mock_constraints['EpGBM']['err_unnormed']
        hist = pop.mock_constraints['EpGBM']['hist_unnormed']
        # Turn into density
        bins_log = np.log10(bins)
        delta_bin = bins_log[1:]-bins_log[:-1]
        N_tot = np.sum(hist)
        hist = hist / (N_tot*delta_bin)
        err = err / (N_tot*delta_bin)

        # Make it look better for plots
        x_plottable, y_mod_plottable = plottable_hist(bins_log[:-1], hist, last_bin_edge=bins_log[-1])
        y_mod_err_plottable = double_array(err)

        # Plot them
        art, = ax.plot(x_plottable, y_mod_plottable, **kwargs)
        ax.fill_between(x_plottable,
                        y_mod_plottable+y_mod_err_plottable,
                        y_mod_plottable-y_mod_err_plottable,
                        alpha=0.6, color=plt.getp(art,'color'), zorder=plt.getp(art,'zorder'))
    ax.set_ylim(0)
    ax.set_xlim(0.5,4)
    ax.set_xlabel(r'log $E_{p {\rm obs}}~\rm{[keV]}$')
    ax.set_ylabel(r'Number density')

    return


def plot_redshift_constraint(ax, pop=None, plot_obs=True, **kwargs):

    if plot_obs:
        bins_eBAT6, obs_lin, err_lin = io.read_constraint(name='eBAT6', density=True, last_bin_edge=6)
        x, x_errp, x_errm = xerr_from_bins(bins_eBAT6)
        plot_obs_property(root_dir/'catalogs/BAT6_cat/eBAT6_cat.txt', header=3,
                          key='redshift',
                          log=False, verbose=True, kde=True,
                          bins=bins_eBAT6,
                          ax=ax, density=True, label=None, alpha=0.8, zorder=0,
                          color='lightgray', errors='coerce')
        ax.errorbar(x, obs_lin, xerr=[x_errm, x_errp], yerr=err_lin, color='k', fmt='none')
        ax.get_legend().remove()

    if pop is not None:
        # Read the mock constraint from the population
        bins = pop.mock_constraints['eBAT6']['bins']
        err = pop.mock_constraints['eBAT6']['err_unnormed']
        hist = pop.mock_constraints['eBAT6']['hist_unnormed']
        # Turn into density
        delta_bin = bins[1:]-bins[:-1]
        N_tot = np.sum(hist)
        hist = hist / (N_tot*delta_bin)
        err = err / (N_tot*delta_bin)

        # Make it look better for plots
        x_plottable, y_mod_plottable = plottable_hist(bins[:-1], hist, last_bin_edge=bins[-1])
        y_mod_err_plottable = double_array(err)

        # Plot them
        art, = ax.plot(x_plottable, y_mod_plottable, **kwargs)
        ax.fill_between(x_plottable,
                        y_mod_plottable+y_mod_err_plottable,
                        y_mod_plottable-y_mod_err_plottable,
                        alpha=0.6, color=plt.getp(art,'color'), zorder=plt.getp(art,'zorder'))
    ax.set_ylim(ymin=0)
    ax.set_xlim(0,6)
    ax.set_xlabel('Redshift(z)')
    ax.set_ylabel(r'Number density')
    return


def plot_SHOALS_z_distr_with_bounds(ax):
    """
        Plot the observed SHOALS cumulative redshift distribution
        while calculating the bounds at the 95% confidence on a given ax
    """
    fname = root_dir/'catalogs/SHOALS_cat/SHOALS_cat.txt'
    z_SHOALS = io.read_data(fname, 2, err=False)
    mask = np.isfinite(z_SHOALS[0])
    z_SHOALS = msc.mask_ndarray(z_SHOALS, mask)
    plot_ndarray_cdf_lim_only(ax, z_SHOALS, arrow_size=0.015)
    st.compute_CDF_bounds_by_MC(sample=z_SHOALS[0],
                                sample_errp=z_SHOALS[1],
                                sample_errm=z_SHOALS[2],
                                sample_ll=z_SHOALS[4],
                                sample_ul=z_SHOALS[3],
                                ul_min_val=np.zeros(len(z_SHOALS[3])),
                                confidence=95.0,
                                positive=True,
                                precision=1000,
                                N_MC=10000,
                                bootstrap=True,
                                verbose=True,
                                show_plot=True,
                                ax=ax,
                                color='k',
                                show_median=True, zorder=0, label='SHOALS observed', linewidth=2)
    ax.set_xlabel('Redshift (z)')
    ax.set_xlim(0,8)
    ax.set_ylim(0,1)
    ax.legend(loc='center right')
    return


def plot_SHOALS_distr(pop, key, SHOALS_file=None, cumul=False, ax=None, plot_obs=True,
    mod_color='C0', mod_label='Model', bins=None, log=False, **kwargs):
    """
        A general function to plot one of two possible distributions for the SHOALS sample:
        Redshift distribution (key='z') or Fluence distribution (key='S_BAT')
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))

    if key == 'z':
        key_mod = 'z'
    elif key == 'S_BAT':
        key_mod = 'pht_flnc_BAT'
    else:
        raise ValueError
    df_obs = io.read_SHOALS_file(SHOALS_file)
    data_obs = df_obs[key].to_numpy()

    try:
        cond = pop.properties['pdet_SHOALS'] == 1
    except KeyError:
        cond = pop.properties['erg_flnc_BAT'] >= 1e-6
    mod = pop.properties[cond][key_mod].to_numpy()

    if log:
        data_obs = np.log10(data_obs)
        mod = np.log10(mod)

    if cumul:
        if plot_obs:
            x_obs, cdf_obs = st.unbinned_empirical_cdf(data_obs)
            ax.plot(x_obs, cdf_obs, color='C5', drawstyle='steps-post', label='SHOALS observed', lw=3)

        x_mod, cdf_mod = st.unbinned_empirical_cdf(mod)
        ax.plot(x_mod, cdf_mod, color=mod_color, drawstyle='steps-post', label=mod_label, **kwargs)

        ax.legend(loc='lower right')
        ax.set_ylabel('CDF')
        ax.set_ylim(0,1)
    else:
        if plot_obs:
            sns.kdeplot(df_obs[key], ax=ax, label='SHOALS observed', color='C5')
            ax.hist(df_obs[key], color='C5', bins=bins, density=True, alpha=0.7)

        sns.kdeplot(mod, ax=ax, label='Model', color=mod_color)
        ax.hist(mod, color=mod_color, bins=bins, density=True, alpha=0.7, **kwargs)

        ax.set_ylabel('PDF')
    return


def plot_Pescalli_2016_GRB_rate(ax, nGRB0=1.3e-9, **kwargs):
    """
        Plot the redshift distribution of LGRBs derived by Pescalli+16
        (https://ui.adsabs.harvard.edu/abs/2016A%26A...587A..40P/abstract)
        DOI : 10.1051/0004-6361/201526760
        The distribution is normalized to its maximum, to convert it to
        physical units use nGRB0 (default value is adjusted to redshift
        distribution of Wanderman & Piran 2010 at z=0).
    """
    filename = root_dir/'catalogs/BAT6_cat/BAT6ext_GRB_formation_rate.txt'

    z = io.read_column(filename, 0) - 1.  # because the original is as 1+z
    distr = io.read_column(filename, 1)
    distr_err = io.read_column(filename, 2)
    norm = nGRB0/distr[0] * 1.6
    # the additional factor comes from the fact that the first value is not
    # exactly as z=0, so we need to correct for that
    ax.errorbar(z, norm*distr, yerr=norm*distr_err, **kwargs)
    ax.fill_between(z, norm*(distr+distr_err), norm*(distr-distr_err), alpha=0.3, color=kwargs['color'])
    return


def plot_constraint(constraint, mod=None, plot_mod=True, plot_obs=True, ax=None):
    """
        Convenience function to quickly plot the various constraints
        used in the population model
    """
    if constraint not in ['Stern', 'EpGBM', 'eBAT6']:
        raise ValueError("constraint must be one of {}".format(['Stern', 'EpGBM', 'eBAT6']))

    if ax is None:
        fig, ax = plt.subplots()

    if not plot_mod:
        mod = None

    if constraint == 'Stern':
        plot_intensity_constraint(ax, pop=mod, plot_obs=plot_obs)
    elif constraint == 'EpGBM':
        plot_spectral_constraint(ax, pop=mod, plot_obs=plot_obs)
    elif constraint == 'eBAT6':
        plot_redshift_constraint(ax, pop=mod, plot_obs=plot_obs)

    ax.legend()
    return


def scatter_incomplete_ndarray(ax, x, y, colormap=None, x_is_log=False, y_is_log=False,
    xlimsize=None, ylimsize=None, capsize=0.0, errorcolor='k', errlw=1.2, z_order=5,
    alpha_errorbar=0.8, edgecolor='k', linewidth=0.8, **kwargs):
    """
    Helper function to easily plot two variables with incomplete data sets (uses masks) and errors, or limits.
    Assumes the data is in the form of the output of read_data (ndarray).
    i.e. : x[0] = data (float)
           x[1] = error plus (float)
           x[2] = error minus (float)
           x[3] = upper limit (bool)
           x[4] = lower limit (bool)
    Colormap input is expected as such:
        Colormap[0] is an array of values assumed to be of same length as x and y.
        Colormap[1] is one value (lower value for colorbar scale)
        Colormap[2] is one value (upper value for colorbar scale)
        Colormap[3] is the argument to put in cmap (ex: 'YlOrRd_r')
    ax : axes object on which to plot.
    Returns scatter plot artist (used to create colorbar afterward).
    """

    # Create masks
    x_mask = np.isfinite(x[0])
    y_mask = np.isfinite(y[0])

    # filters data
    x_to_plot = x[0][x_mask & y_mask]
    xerr = np.asarray([np.nan_to_num(x[2][x_mask & y_mask]), np.nan_to_num(x[1][x_mask & y_mask])])
    xuplims = x[3][x_mask & y_mask]
    xlolims = x[4][x_mask & y_mask]
    xmin = np.min(x_to_plot-xerr[0])
    xmax = np.max(x_to_plot+xerr[0])

    if xlimsize is None:
        xlimsize = (xmax-xmin)/10.

    for i in range(len(x_to_plot)):
        if x_is_log:
            if xuplims[i]:
                xerr[0][i] = 0.5 * x_to_plot[i]  # lower error for upper limit (arrow pointing down)
            if xlolims[i]:
                xerr[1][i] = 1.5 * x_to_plot[i]  # upper error for lower limit (arrow pointing up)
            # The coefficients were chosen to give sizeable arrows in log scale
        else:
            if xuplims[i]:
                xerr[0][i] = xlimsize  # lower error for upper limit (arrow pointing down)
            if xlolims[i]:
                xerr[1][i] = xlimsize  # upper error for lower limit (arrow pointing up)

    y_to_plot = y[0][x_mask & y_mask]
    yerr = np.asarray([np.nan_to_num(y[2][x_mask & y_mask]), np.nan_to_num(y[1][x_mask & y_mask])])
    yuplims = y[3][x_mask & y_mask]
    ylolims = y[4][x_mask & y_mask]
    ymin = np.min(y_to_plot-yerr[0])
    ymax = np.max(y_to_plot+yerr[0])
    if ylimsize is None:
        ylimsize = (ymax-ymin)/10.

    for i in range(len(y_to_plot)):
        if y_is_log:
            if yuplims[i]:
                yerr[0][i] = 0.5 * y_to_plot[i]  # lower error for upper limit (arrow pointing down)
            if ylolims[i]:
                yerr[1][i] = 1.5 * y_to_plot[i]  # upper error for lower limit (arrow pointing up)
            # The coefficients were chosen to give sizeable arrows in log scale
        else:
            if yuplims[i]:
                yerr[0][i] = ylimsize  # lower error for upper limit (arrow pointing down)
            if ylolims[i]:
                yerr[1][i] = ylimsize  # upper error for lower limit (arrow pointing up)

    # Plotting
    ax.errorbar(x_to_plot, y_to_plot, xerr=xerr, yerr=yerr,
                xuplims=xuplims, xlolims=xlolims, uplims=yuplims, lolims=ylolims,
                capsize=capsize, color=errorcolor, marker=None, linewidth=errlw, fmt='.',
                alpha=alpha_errorbar, zorder=z_order)
    if colormap is not None:
        norm = Normalize(vmin=colormap[1], vmax=colormap[2])  # limits to the colorbar if colormap is used
        scatterplot = ax.scatter(x_to_plot, y_to_plot, c=colormap[0][x_mask & y_mask],
                                 cmap=colormap[3], norm=norm, zorder=z_order+1,
                                 edgecolor=edgecolor, linewidth=linewidth, **kwargs)
    else:
        scatterplot = ax.scatter(x_to_plot, y_to_plot, zorder=z_order+1, edgecolor=edgecolor,
                                 linewidth=linewidth, **kwargs)

    return scatterplot


def fig_marg(figsize=(10, 8), cb=True, **kwargs):
    """
        Create a figure with marginal plots on the top and the left.
        Optional: cb stands for colorbar
    """

    fig = plt.figure(figsize=figsize, tight_layout=False, **kwargs)

    if cb:
        ax_center = fig.add_axes([0.25, 0.10, 0.60, 0.70])
        ax_left = fig.add_axes([0.10, 0.10, 0.15, 0.70], sharey=ax_center)
        ax_top = fig.add_axes([0.25, 0.80, 0.60, 0.15], sharex=ax_center)
        ax_cb = fig.add_axes([0.86, 0.1, 0.03, 0.70])
    else:
        ax_center = fig.add_axes([0.25, 0.10, 0.65, 0.70])
        ax_left = fig.add_axes([0.10, 0.10, 0.15, 0.70], sharey=ax_center)
        ax_top = fig.add_axes([0.25, 0.80, 0.65, 0.15], sharex=ax_center)

    ax_left.invert_xaxis()
    ax_center.tick_params(axis='y', which='both', labelleft=False)
    ax_left.tick_params(axis='x', which='both', labelbottom=False)
    ax_top.tick_params(axis='both', which='both', labelleft=False, labelbottom=False)
    axes = {'center':ax_center,
            'left':ax_left,
            'top':ax_top}
    if cb:
        ax_cb.yaxis.set_label_position('right')
        ax_cb.yaxis.set_ticks_position('right')
        axes['cb'] = ax_cb
    return fig, axes


def cool_hist2d(x, y, c=None, mode='scatter', cb=True, fig=None, figsize=(10, 8),
    xlabel=None, ylabel=None, cblabel=None,
    plot_left_kdeplot=True, plot_top_kdeplot=True, plot_left_hist=True, plot_top_hist=True,
    left_kdeplot_kwargs={'color': 'k', 'label': 'KDE'},
    top_kdeplot_kwargs={'color': 'k', 'label': 'KDE'},
    left_hist_kwargs={'label': None, 'bins': 20, 'color': 'lightgrey', 'edgecolor':'k'},
    top_hist_kwargs={'label': None, 'bins': 20, 'color': 'lightgrey', 'edgecolor':'k'},
    cbar_kwargs={}, hist2d_kwargs={}, kde2d_kwargs={}, **kwargs):
    """
        Convenience function to plot a 2D x y plane.
        There are 3 values possible for 'mode':
            - scatter: scatter plot
            - hist2d: a 2D histogram of the data
            - kde2D: a 2D Gaussian Kernel Density Estimation of the data
        c is an optional argument which can be used to color the points
        in the plane with a third property in the 'scatter' mode.
    """
    if fig is None:
        fig, axes = fig_marg(figsize=figsize, cb=cb)
    else:
        ax_list = fig.axes
        axes = {}
        axes['center'] = ax_list[0]
        axes['left'] = ax_list[1]
        axes['top'] = ax_list[2]
        if cb:
            axes['cb'] = ax_list[3]

    if mode == 'scatter':
        art = axes['center'].scatter(x, y, c=c, **kwargs)
    elif mode == 'hist2d':
        corner.hist2d(x, y, ax=axes['center'], **hist2d_kwargs)
    elif mode == 'kde2d':
        sns.kdeplot(x, y, ax=axes['center'], **kde2d_kwargs)
    else:
        raise ValueError('Wrong input for mode in cool_hist2d. Please chose "scatter" or "hist2d"')

    if plot_left_kdeplot:
        sns.kdeplot(y, ax=axes['left'], vertical=True, **left_kdeplot_kwargs)
    if plot_top_kdeplot:
        sns.kdeplot(x, ax=axes['top'], **top_kdeplot_kwargs)
    if plot_left_hist:
        axes['left'].hist(y, orientation='horizontal', density=True, **left_hist_kwargs)
    if plot_top_hist:
        axes['top'].hist(x, density=True, **top_hist_kwargs)

    axes['center'].set_xlabel(xlabel)
    axes['left'].set_ylabel(ylabel)

    if cb and c is not None:
        cb = fig.colorbar(art, cax=axes['cb'], **cbar_kwargs)
    if cb and cblabel is not None:
        cb.set_label(cblabel)
    return fig, axes


def plot_obs_property(fname, key, func=None, func_args={}, ax=None, log=False, kde=True,
    header=2, verbose=False, debug=False, errors='raise', **kwargs):
    """
        Convenience function to quickly plot an observed sample from a given file name.
        A function to filter or cut the sample can be passed as func.
    """
    # Read the entire file
    df_obs = pd.read_csv(fname, sep='|', header=header, low_memory=False)
    # Strip the colum names to remove whitespaces
    df_obs.rename(columns=lambda x:x.strip(), inplace=True)
    # Activate debug to check the column names
    if debug:
        for i,col in enumerate(df_obs.columns):
            print(i,col)

    # Apply function to the data
    if func is None:
        df_prop = pd.to_numeric(df_obs[key], errors=errors)
    # If func is a list, iterate through the list and apply each function
    elif isinstance(func, list):
        if not isinstance(func_args, list):
            raise ValueError
        df_prop = df_obs.copy()
        for i, func_i in enumerate(func):
            df_prop = func_i(df_prop, **func_args[i])
        df_prop = pd.to_numeric(df_prop[key], errors=errors)
    else:
        df_prop = func(df_obs.copy(), **func_args)
        df_prop = pd.to_numeric(df_prop[key], errors=errors)

    if log:
        df_prop = pd.to_numeric(df_prop, errors=errors)
        df_prop = np.log10(df_prop)
    if verbose:
        print("Sample size :{}".format(len(df_prop.dropna())))
    if ax is None:
        ax = plt.gca()

    ax.hist(df_prop, **kwargs)
    if kde:
        sns.kdeplot(df_prop.dropna(), ax=ax, color='k', linewidth=2, label='KDE')
    return df_prop


def plot_logNlogP(fname, key, func=None, func_args={}, bins=None, sep=None, hdr=None, density=False,
                  cols=None, ax=None, eff=False, errors='raise', y_log=False, bins_log=True, factor=None,
                  Stern_first_bin_corr=False, cat_duration=False, **kwargs):
    """
        Plot the logN logP diagram for a given catalog
    """

    # Read data
    df_obs = pd.read_csv(fname, sep=sep, header=hdr, names=cols, low_memory=False)
    # Strip the colum names to remove whitespaces
    df_obs.rename(columns=lambda x:x.strip(), inplace=True)
    print("Freshly read N_data: {}".format(len(df_obs.index)))
    if cat_duration:
        msc.calc_cat_duration(fname, verbose=True)
    # Apply function to the data
    if func is None:
        df_prop = pd.to_numeric(df_obs[key], errors=errors)
    # If func is a list, iterate through the list and apply each function
    elif isinstance(func, list):
        if not isinstance(func_args, list):
            raise ValueError
        df_prop = df_obs.copy()
        for i, func_i in enumerate(func):
            df_prop = func_i(df_prop, **func_args[i])
        df_prop = pd.to_numeric(df_prop[key], errors=errors)
    else:
        df_prop = func(df_obs.copy(), **func_args)
        df_prop = pd.to_numeric(df_prop[key], errors=errors)

    print("After filtering, N_data (excluding NaNs): {}".format(df_prop.count()))
    if y_log:
        df_prop = np.log10(df_prop)

    if bins is None:
        bins, _u, _u = io.read_logRlogN()
    x, x_errp, x_errm = xerr_from_bins(bins, logscale=bins_log)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))
    hist, bins = np.histogram(df_prop, bins=bins)
    hist_err = np.sqrt(hist)

    if eff:
        hist = hist / msc.efficiency_correction_Stern(x)
        hist_err = hist_err / msc.efficiency_correction_Stern(x)
    if density:
        N_tot = np.sum(hist)
    else:
        N_tot = 1
    delta_bins = np.log10(bins[1:]/bins[:-1])
    hist = hist/(N_tot*delta_bins)
    hist_err = hist_err/(N_tot*delta_bins)

    if factor is not None:
        hist = hist/factor
        hist_err = hist_err/factor
    if Stern_first_bin_corr:
        hist[0], hist_err[0], _ = msc.log_to_lin(3.069180, 9.691000e-02)
    ax.errorbar(x, hist, xerr=[x_errm, x_errp], yerr=hist_err,
                fmt='none', **kwargs)

    return hist, hist_err


def plot_eBAT6_EpL(fname=None, axes=None, fname_lim=None, mini_cax=False, kde=False,
    show_relation=True):
    if fname is None:
        fname = root_dir/'catalogs/BAT6_cat/eBAT6_cat.txt'
    if fname_lim is None:
        fname_lim = root_dir/'resources/BAT6_detection_limit_EpL_plane.txt'

    if axes is None:
        fig, axes = fig_marg(figsize=(10, 8), cb=True)
    else:
        fig = plt.gcf()

    if mini_cax:
        m_cax = fig.add_axes([0.780, 0.125, 0.015, 0.20])
        axes['mini_cb'] = m_cax
        fsize = 12
    else:
        axes['mini_cb'] = axes['cb']
        fsize = 25

    # Retrieve observational data
    Ep_obs = io.read_data(fname, 13, stripper='|', splitter='|', single_err=True)
    L_obs = io.read_data(fname, 15, stripper='|', splitter='|', single_err=True)
    z_obs = io.read_data(fname, 1, stripper='|', splitter='|', err=False)
    L_obs[0:3] = L_obs[0:3] * 10**51
    mask = np.isfinite(L_obs[0]) & np.isfinite(Ep_obs[0])
    z_obs_masked = msc.mask_ndarray(z_obs, mask)
    Ep_obs_masked = msc.mask_ndarray(Ep_obs, mask)
    L_obs_masked = msc.mask_ndarray(L_obs, mask)
    Ep_obs_masked = msc.lin_to_log_ndarray(Ep_obs_masked)
    L_obs_masked = msc.lin_to_log_ndarray(L_obs_masked)

    # Plot them
    art = scatter_incomplete_ndarray(axes['center'], L_obs_masked, Ep_obs_masked,
                                     colormap=[z_obs_masked[0],0.,6.,'YlOrRd'],
                                     marker='o', s=100, edgecolor='k', linewidth=0.8,
                                     label='eBAT6 observed data')
    if kde:
        sns.kdeplot(L_obs_masked[0], Ep_obs_masked[0], ax=axes['center'], n_levels=5,
                    shade=True, shade_lowest=False, cmap='Greys', alpha=0.8)
    cb2 = fig.colorbar(art, cax=axes['mini_cb'])
    cb2.set_label('Redshift (z)', **{'size':fsize})
    cb2.ax.tick_params(labelsize=20)
    cb2.ax.tick_params(which='minor', right=False)

    # Detection limit of BAT6
    Ep_lim = io.read_column(fname_lim, 0, splitter='\t')
    z_lim = np.asarray([0.3, 2., 5.])
    L_lim = np.zeros((len(z_lim), len(Ep_lim)))
    ls = ['--', '-.',':']
    for i in range(len(z_lim)):
        L_lim[i] = io.read_column(fname_lim, i+1, splitter='\t')
        axes['center'].plot(np.log10(L_lim[i]), np.log10(Ep_lim), ls=ls[i], c='k', lw=2)

    # Relation of Pescalli+16
    if show_relation:
        logL_Pesc = np.linspace(49, 55)
        logEp_Pesc = -25.6 + 0.54*logL_Pesc
        sigma_Pesc = 0.28 * np.sqrt(1.+0.54**2)  # to adjust for perpendicular scatter
        # Ep0 = 309. keV, alpha_amati = 0.54, sigma = 0.28
        art, = axes['center'].plot(logL_Pesc, logEp_Pesc, lw=2.5, color='C0')
        color_Pesc = plt.getp(art,'color')
        axes['center'].plot(logL_Pesc, logEp_Pesc+sigma_Pesc, ls='--', lw=1.5, color=color_Pesc)
        axes['center'].plot(logL_Pesc, logEp_Pesc-sigma_Pesc, ls='--', lw=1.5, color=color_Pesc)

    # Left histogram
    bins_Ep = np.linspace(1,4,31)
    axes['left'].hist(Ep_obs_masked[0], bins=bins_Ep, orientation='horizontal', label='eBAT6 observed',
                      edgecolor='k', linewidth=0.5, density=True, color='lightgray', alpha=0.8)
    sns.kdeplot(Ep_obs_masked[0], ax=axes['left'], vertical=True, color='k', linewidth=2)
    axes['left'].invert_xaxis()
    axes['left'].autoscale(True, axis='x')
    # Top histogram
    bins_L = np.linspace(49,55,31)
    axes['top'].hist(L_obs_masked[0], bins=bins_L, label='eBAT6 observed',
                     edgecolor='k',linewidth=0.5, density=True, color='lightgray', alpha=0.8)
    sns.kdeplot(L_obs_masked[0], ax=axes['top'], label='KDE', color='k', linewidth=2)
    axes['top'].get_legend().remove()
    axes['top'].autoscale(True, axis='y')

    # Axes formatting
    axes['left'].set_ylabel(r'log Peak Energy $\rm{[keV]}$')
    axes['center'].set_xlabel(r'log Luminosity $\rm{[erg/s]}$')
    axes['center'].set_xlim(49.3,54.8)
    axes['center'].set_ylim(1.3, 3.7)
    axes['center'].tick_params(labelsize=23)
    axes['left'].tick_params(labelsize=23)

    return


def plot_CDF_with_bounds(bins_mid, median, lower, upper, ax, **kwargs):
    """
        Convenience function to elegantly plot a cumulative distribution
        function and the confidence bounds around it.
    """

    artist, = ax.plot(bins_mid, median, drawstyle='steps-mid', **kwargs)

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
    ax.add_line(bottom_line)
    top_line = Line2D([xmax, xmax+(xmax-xmin)], [median[imax],median[imax]],
                      color=plot_color, linestyle=plot_ls, linewidth=plot_lw, zorder=plot_zorder)
    ax.add_line(top_line)
    ax.fill_between(bins_mid, lower, upper, step='mid', color=plot_color, alpha=0.3, zorder=plot_zorder-1)
    ax.plot(bins_mid, lower, drawstyle='steps-mid',lw=0.7, c=plot_color, zorder=plot_zorder-1)
    ax.plot(bins_mid, upper, drawstyle='steps-mid',lw=0.7, c=plot_color, zorder=plot_zorder-1)
    ax.legend()
    return


def plot_CDFs_and_KS_results(bins, med1, med2, lw1, lw2, up1, up2, D_stat, p_value, confidence,
    pfrac=None, label1=None, label2=None, color1=None, color2=None):
    fig = plt.figure(figsize=(10,8), tight_layout=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1], width_ratios=[1,1])
    ax = fig.add_subplot(gs[0,:])
    axD = fig.add_subplot(gs[1,0])
    axp = fig.add_subplot(gs[1,1])

    ax.set_ylabel('CDF')
    axD.set_xlabel('D-stat')
    axp.set_xlabel('log(p-value)')
    axD.hist(D_stat, bins=20, color='lightgray', density=True)
    if pfrac is None:
        pfrac = len(p_value[np.where(p_value < (1.-confidence/100.))[0]])/len(p_value)

    plabel = r"{:.2f} \% below {:.2f}".format(100*pfrac, 1.-confidence/100.)
    axp.hist(np.log10(p_value), bins=20, color='lightgray', density=True)
    axp.axvline(np.log10(1.-confidence/100.), ls='--', color='k', label=plabel)
    axp.legend()

    bins_mid = 0.5*(bins[1:]+bins[:-1])
    plot_CDF_with_bounds(bins_mid, med1, lw1, up1, ax=ax, label=label1, color=color1)
    plot_CDF_with_bounds(bins_mid, med2, lw2, up2, ax=ax, label=label2, color=color2)
    ax.set_ylim(0,1)
    return


def plot_ndarray_unbinned_cdf_with_limits(ax, data, x_is_log=False, arrow_size=None, **kwargs):
    """
        Helper function to quickly plot a cumulative distribution including limits.
        The input is expected to come from read_data
        i.e. : x[0] = data (float)
           x[1] = error plus (float)
           x[2] = error minus (float)
           x[3] = upper limit (bool)
           x[4] = lower limit (bool)
        Returns the artist associated with the plotted line
    """
    mask = np.isfinite(data[0])
    masked_data = msc.mask_ndarray(data, mask)
    sorted_data = msc.sort_ndarray(masked_data)
    w_ul = np.where(sorted_data[3] == 1)[0]
    w_ll = np.where(sorted_data[4] == 1)[0]

    _unused, ECDF = st.unbinned_empirical_cdf(sorted_data[0])

    artist, = ax.plot(sorted_data[0], ECDF, drawstyle='steps-post', **kwargs)

    # Add the first and last line to make the plot look better
    plot_color = plt.getp(artist, 'color')
    plot_ls = plt.getp(artist, 'linestyle')
    plot_lw = plt.getp(artist, 'linewidth')
    xmin = sorted_data[0,0]
    xmax = sorted_data[0,-1]
    bottom_line = Line2D([xmin,xmin], [0,ECDF[0]],
                         color=plot_color, linestyle=plot_ls, linewidth=plot_lw)
    ax.add_line(bottom_line)
    top_line = Line2D([xmax, (2*xmax-sorted_data[0,-2])], [ECDF[-1],ECDF[-1]],
                      color=plot_color, linestyle=plot_ls, linewidth=plot_lw)
    ax.add_line(top_line)

    if x_is_log:
        if arrow_size is None:
            arrow_size = [0.3, 0.8]
        arrow_length = np.ones(len(sorted_data[0]))
        arrow_length[w_ul] = sorted_data[0, w_ul] * arrow_size[0]
        arrow_length[w_ll] = sorted_data[0, w_ll] * arrow_size[1]
    else:
        if arrow_size is None:
            arrow_size = 0.1
        arrow_length = np.ones(len(sorted_data[0]))*(xmax-xmin) * arrow_size
    # Add upper limits
    for i in range(len(w_ul)):
        if w_ul[i] > 0:
            ax.arrow(sorted_data[0, w_ul[i]], 0.5*(ECDF[w_ul[i]]+ECDF[w_ul[i]-1]), -arrow_length[w_ul[i]], 0,
                     fc=plot_color, ec=plot_color, head_width=0.01, head_length=0.1*arrow_length[w_ul[i]])
        else:
            ax.arrow(sorted_data[0, w_ul[i]], 0.5*(ECDF[w_ul[i]]), -arrow_length[w_ul[i]], 0,
                     fc=plot_color, ec=plot_color, head_width=0.01, head_length=0.1*arrow_length[w_ul[i]])
    # and lower limits
    for i in range(len(w_ll)):
        if w_ll[i] > 0:
            ax.arrow(sorted_data[0, w_ll[i]], 0.5*(ECDF[w_ll[i]]+ECDF[w_ll[i]-1]), arrow_length[w_ll[i]], 0,
                     fc=plot_color, ec=plot_color, head_width=0.01, head_length=0.1*arrow_length[w_ll[i]])
        else:
            ax.arrow(sorted_data[0, w_ll[i]], 0.5*(ECDF[w_ll[i]]), arrow_length[w_ll[i]], 0,
                     fc=plot_color, ec=plot_color, head_width=0.01, head_length=0.1*arrow_length[w_ll[i]])

    return artist


def plot_ndarray_cdf_lim_only(ax, ndarr, x_is_log=False, arrow_size=None, color='k', alpha=0.8,
    loc='bottom', **kwargs):

    xmin = ndarr[0].min()
    xmax = ndarr[0].max()
    w_ul = np.where(ndarr[3] == 1)[0]
    w_ll = np.where(ndarr[4] == 1)[0]
    if x_is_log:
        if arrow_size is None:
            arrow_size = [0.3, 0.8]
        arrow_length = np.ones(len(ndarr[0]))
        arrow_length[w_ul] = ndarr[0, w_ul] * arrow_size[0]
        arrow_length[w_ll] = ndarr[0, w_ll] * arrow_size[1]
    else:
        if arrow_size is None:
            arrow_size = 0.1
        arrow_length = np.ones(len(ndarr[0]))*(xmax-xmin) * arrow_size
    if loc == 'bottom':
        y_tail_start = 0.01
        y_tail_stop = 0.03
    elif loc == 'top':
        y_tail_start = 0.97
        y_tail_stop = 0.99

    for i in range(len(w_ul)):
        # upper limit
        l2D = ax.axvline(ndarr[0,w_ul[i]], y_tail_start, y_tail_stop, color=color, linewidth=2)
        ax.arrow(ndarr[0, w_ul[i]], np.average(l2D.get_ydata()), -arrow_length[w_ul[i]], 0,
                 fc=color, ec=color, head_width=0.01, head_length=0.4*arrow_length[w_ul[i]], alpha=alpha, **kwargs)
    for i in range(len(w_ll)):
        l2D = ax.axvline(ndarr[0,w_ll[i]], y_tail_start, y_tail_stop, color=color, linewidth=2)
        ax.arrow(ndarr[0, w_ll[i]], np.average(l2D.get_ydata()), arrow_length[w_ll[i]], 0,
                 fc=color, ec=color, head_width=0.01, head_length=0.4*arrow_length[w_ll[i]], alpha=alpha, **kwargs)

    return


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


def xerr_from_bins(bins, logscale=False):
    """
        Convenience function to easily calculate the x_err from bins to
        represent a histogram as crosses on a plot
        returns:
        - x the middle point of the bin in the x direction
        - x_errp the extent of the bin from x up to the right bin edge
        - x_errm the extent of the bin from x down to the left bin edge
    """
    if logscale:
        x = np.sqrt(bins[1:] * bins[:-1])
    else:
        x = 0.5 * (bins[1:] + bins[:-1])

    x_errp = bins[1:] - x
    x_errm = x - bins[:-1]

    return x, x_errp, x_errm
