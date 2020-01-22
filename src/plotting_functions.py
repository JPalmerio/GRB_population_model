import corner
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import observational_constraints as obs
import stats as st
from io_grb_pop import root_dir
import io_grb_pop as io
from matplotlib import gridspec

try:
    plt.style.use('presentation')
except Exception:
    pass


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


def plot_EpGBM(fname=None, log=True, density=False, ax=None, **kwargs):
    """
        Convenience function to plot the EpGBM constraint on a given ax
    """
    if fname is None:
        fname = root_dir/'observational_constraints/Ep_GBM_for_plotting.txt'
    if ax is None:
        ax = plt.gca()
    bins, obs_lin, err_lin = obs.create_EpGBM_hist(fname, density=density)

    x, x_errp, x_errm = xerr_from_bins(bins, log=log)

    ax.errorbar(x, obs_lin, xerr=[x_errm, x_errp], yerr=err_lin, **kwargs)

    return


def plot_SHOALS_distr(df_mod, key, SHOALS_file=None, cumul=False, fig=None, plot_obs=True,
                      mod_color='C0', mod_label='Model', bins=None):

    if fig is None:
        fig, ax = plt.subplots(figsize=(10,8))
    else:
        ax = fig.axes[0]

    df_obs = io.read_SHOALS_file(SHOALS_file)
    if cumul:
        if plot_obs:
            if key == 'S_BAT':
                data_obs = np.log10(df_obs[key].to_numpy())
            else:
                data_obs = df_obs[key].to_numpy()

            x_obs, cdf_obs = st.unbinned_empirical_cdf(data_obs)
            ax.plot(x_obs, cdf_obs, color='C5', drawstyle='steps-post', label='SHOALS observed', lw=3)

        x_mod, cdf_mod = st.unbinned_empirical_cdf(df_mod.to_numpy())
        ax.plot(x_mod, cdf_mod, color=mod_color, drawstyle='steps-post', label=mod_label)

        ax.legend()
        ax.set_ylabel('CDF')
        ax.set_ylim(0,1)
    else:
        if (bins is None) and (key == 'z'):
            bins = np.linspace(0,7,15)
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


def plot_constraint(constraint, bins, mod, obs, err, plot_mod=True, plot_obs=True, log=True, ax=None):
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


def fig_marg(figsize=(7, 6), cb=True):
    fig = plt.figure(figsize=figsize)

    if cb:
        ax_center = fig.add_axes([0.25, 0.10, 0.60, 0.70])
        ax_left = fig.add_axes([0.10, 0.10, 0.15, 0.70], sharey=ax_center)
        ax_top = fig.add_axes([0.25, 0.80, 0.60, 0.15], sharex=ax_center)
        ax_cb = fig.add_axes([0.86, 0.1, 0.03, 0.75])
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


def cool_hist2d(x, y, c=None, mode='scatter', xlabel=None, ylabel=None, cb=True, cblabel=None, fig=None, figsize=(9, 7),
                plot_left_kdeplot=True, plot_top_kdeplot=True,
                plot_left_hist=True, plot_top_hist=True,
                left_kdeplot_kwargs={'color': 'k', 'label': 'KDE'},
                top_kdeplot_kwargs={'color': 'k', 'label': 'KDE'},
                left_hist_kwargs={'label': None,
                                  'bins': 20,
                                  'color': 'lightgrey',
                                  'edgecolor':'k'},
                top_hist_kwargs={'label': None,
                                 'bins': 20,
                                 'color': 'lightgrey',
                                 'edgecolor':'k'},
                cbar_kwargs={}, hist2d_kwargs={},
                **kwargs):
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
        art = axes['center'].scatter(x, y, c=c, edgecolor='k', **kwargs)
    elif mode == 'hist2d':
        corner.hist2d(x, y, ax=axes['center'], **hist2d_kwargs)
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
    bottom_line = matplotlib.lines.Line2D([xmin,xmin], [0,median[imin]],
                                          color=plot_color, linestyle=plot_ls, linewidth=plot_lw, zorder=plot_zorder)
    ax.add_line(bottom_line)
    top_line = matplotlib.lines.Line2D([xmax, xmax+(xmax-xmin)], [median[imax],median[imax]],
                                       color=plot_color, linestyle=plot_ls, linewidth=plot_lw, zorder=plot_zorder)
    ax.add_line(top_line)
    ax.fill_between(bins_mid, lower, upper, step='mid', color=plot_color, alpha=0.3)
    ax.plot(bins_mid, lower, drawstyle='steps-mid',lw=0.7, c=plot_color)
    ax.plot(bins_mid, upper, drawstyle='steps-mid',lw=0.7, c=plot_color)
    ax.legend()
    return


def plot_CDFs_and_KS_results(bins, med1, med2, lw1, lw2, up1, up2, D_stat, p_value, confidence,
                             pfrac=None, label1=None, label2=None):
    fig = plt.figure(figsize=(10,8),tight_layout=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1], width_ratios=[1,1])
    ax = fig.add_subplot(gs[0,:])
    axD = fig.add_subplot(gs[1,0])
    axp = fig.add_subplot(gs[1,1])

    ax.set_ylabel('CDF')
    axD.set_xlabel('D-stat')
    axp.set_xlabel('log(p-value)')
    axD.hist(D_stat, bins=20, color='C2', density=True)
    if pfrac is None:
        pfrac = len(p_value[np.where(p_value < (1.-confidence/100.))[0]])/len(p_value)

    plabel = r"{:.2f} \% below {:.2f}".format(100*pfrac, 1.-confidence/100.)
    axp.hist(np.log10(p_value), bins=20, color='C2', density=True)
    axp.axvline(np.log10(1.-confidence/100.), ls='--', color='gray', label=plabel)
    axp.legend()

    bins_mid = 0.5*(bins[1:]+bins[:-1])
    plot_CDF_with_bounds(bins_mid, med1, lw1, up1, ax=ax, label=label1)
    plot_CDF_with_bounds(bins_mid, med2, lw2, up2, ax=ax, label=label2)
    ax.set_ylim(0,1)
    return
