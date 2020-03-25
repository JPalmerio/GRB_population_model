import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import yaml

filename = '../model_outputs/190329_test_ab/GRB_Properties'
filename_old = '../model_outputs/190329_test_ab_old/GRB_Properties'
with open('../init/samples.yml', 'r') as f:
    samples = yaml.safe_load(f)


def examine_pop(filename, fig_prop=None, fig_pflx=None, fig_samp=None, label=None, **kwargs):
    df = pd.read_msgpack(filename)

    if fig_prop is None:
        fig_prop, axes_prop = plt.subplots(2, 4, figsize=(14, 6), tight_layout=True)
    else:
        axes_prop = np.array(fig_prop.get_axes())

    prop_name = ['z', 'L', 'Ep', 't90', 'Cvar', 'alpha', 'beta', 'ktild']
    are_log = [False, True, True, True, True, False, False, False]
    bins_lims = [(0, 20), (49, 55), (0, 5), (-0.5, 3.5), (-2, 0), (-1, 2), (2, 25), (0, 15)]

    for ax, prop, is_log, bins_lim in zip(axes_prop.flatten(), prop_name, are_log, bins_lims):
        if is_log:
            ax.hist(np.log10(df[prop]), bins=np.linspace(bins_lim[0], bins_lim[1], 20), label=label, density=True, **kwargs)
            ax.set_xlabel(f"log({prop})")
        else:
            ax.hist(df[prop], bins=np.linspace(bins_lim[0], bins_lim[1], 20), label=label, density=True, **kwargs)
            ax.set_xlabel(f"{prop}")
        ax.set_yscale('log')
        ax.legend()

    pflx_name = [key for key in df.columns if (('pflx_' in key) or ('cts_' in key))]
    if fig_pflx is None:
        fig_pflx, axes_pflx = plt.subplots(1, len(pflx_name), figsize=(3.5*len(pflx_name), 4.5), tight_layout=True)
    else:
        axes_pflx = np.array(fig_pflx.get_axes())
    for ax, pflx in zip(axes_pflx.flatten(), pflx_name):
        ax.hist(np.log10(df[pflx]), bins=np.arange(-20, 6, 0.5), label=label, density=True, **kwargs)
        ax.set_xlabel(f"log({pflx})")
        ax.set_xlim(np.log10(df[pflx].min()), np.log10(df[pflx].max()))
        ax.set_yscale('log')
        ax.legend()

    samp_name = [key for key in df.columns if ('pdet_' in key)]
    if fig_samp is None:
        fig_samp, axes_samp = plt.subplots(1, len(samp_name), figsize=(3.5*len(samp_name), 4.5), tight_layout=True)
    else:
        axes_samp = np.array(fig_samp.get_axes())
    for ax, samp in zip(axes_samp.flatten(), samp_name):
        pflx_min = np.log10(float(samples[samp.split('_')[1]]['pflx_min']))
        instrum = '_'.join(['pflx', samples[samp.split('_')[1]]['instrument']])
        ax.set_xlabel(f"log(pflx {'_'.join(samp.split('_')[1:])})")
        if 'ECLAIRs' in samp:
            instrum = '_'.join(['cts', samples[samp.split('_')[1]]['instrument']])
            ax.set_xlabel(f"log(cts {'_'.join(samp.split('_')[1:])})")
        ax.hist(np.log10(df[instrum]), weights=df[samp], bins=np.arange(pflx_min, 7, 0.15), label=label, density=True, **kwargs)
        ax.set_xlim(pflx_min, np.log10(df[instrum].max()))
        ax.set_yscale('log')
        ax.legend()

    return fig_prop, fig_pflx, fig_samp


fig_prop, fig_pflx, fig_samp = examine_pop(filename, label='New', zorder=6, edgecolor='k', linewidth=0.5, alpha=0.8)
fig_prop, fig_pflx, fig_samp = examine_pop(filename_old, fig_prop, fig_pflx, fig_samp, label='Old', zorder=5, edgecolor='k', linewidth=0.5, alpha=0.8)

plt.show()
