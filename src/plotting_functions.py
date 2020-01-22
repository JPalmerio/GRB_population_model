import corner
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec

try:
    plt.style.use('presentation')
except Exception:
    pass


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
