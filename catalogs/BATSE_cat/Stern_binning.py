import sys
import platform
if platform.system() == 'Linux':
    sys.path.insert(0,'/nethome/palmerio/Dropbox/Plotting_GUI/Src')
elif platform.system() == 'Darwin': 
    sys.path.insert(0,'/Users/palmerio/Dropbox/Plotting_GUI/Src')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
import plotting_functions as pf
from pathlib import Path
root_dir = Path('/Users/palmerio/Science_projects/')


def efficiency_correction_Stern(pflx, c_e0=0.097, nu=2.34, norm=0.7):
    """
        The efficiency function of BATSE for detecting GRBs as a function of peak flux, derived by Stern+01
        c_e0 is in [counts/s/cm2]
        pflx is in [ph/s/cm2]
    """
    c_e = pflx * 0.75  # the conversion factor from counts to pflx comes from the Stern+01 paper as well, figure 7.
    return norm * (1.0 - np.exp(-(c_e/c_e0)**2) )**nu


T_min = 2.05   # secondes
P23_min = 0.01  # ph/cm2/s
N_bin_min = 1  # Nb of bins with flx >= 50% of pflx -> more robust than T90>2s because of uncertainties on T90 measurement (see Stern+01)
verbose = 1
file_Stern = str(root_dir/'GRB_population_model/observational_constraints/lognlogp.stern.dat')
text_size = 22
font = {'fontname':'Serif', 'size':text_size}
run = '170704_Amati_zevol'


#################### Stern #######################
filename_S = 'Stern_cat.txt'
name_S = pf.read_column(filename_S, 0, dtype=str, splitter='\t')
T90_S = pf.read_column(filename_S, 7, splitter='\t')
N_bin_S = pf.read_column(filename_S, 8, splitter='\t')
peak_flux_1024_S = pf.read_column(filename_S, 3, splitter='\t')
peak_flux_1024_mask_S = np.isfinite(peak_flux_1024_S)
final_peak_flux_1024_S = peak_flux_1024_S[peak_flux_1024_mask_S]
final_long_peak_flux_1024_S = peak_flux_1024_S[N_bin_S > N_bin_min]

# Stern
bins = pf.read_column(file_Stern, 0, array=False)
hist_Stern = 10**pf.read_column(file_Stern, 1)
hist_Stern_err = 10**pf.read_column(file_Stern, 2)
bins.append(1.8)
bins = 10**np.asarray(bins)/0.75
bins_middle = np.zeros(len(bins)-1)
bins_errp = np.zeros(len(bins_middle))
bins_errm = np.zeros(len(bins_middle))
hist_S, bins_S = np.histogram(final_long_peak_flux_1024_S, bins=bins)
hist_S = np.asarray(hist_S, dtype=float)
hist_S_err = np.sqrt(hist_S)


for i in range(len(bins_middle)):
    bins_middle[i] = np.sqrt(bins[i]*bins[i+1])
    bins_errp[i] = bins[i+1] - bins_middle[i]
    bins_errm[i] = bins_middle[i] - bins[i]
    hist_S[i]         /= np.log10( bins[i+1] / bins[i] )#*0.67*4*np.pi*6.54
    hist_S_err[i]     /= np.log10( bins[i+1] / bins[i] )#*0.67*4*np.pi*6.54
hist_S /= efficiency_correction_Stern(bins_middle)
hist_S_err /= efficiency_correction_Stern(bins_middle)


def decide_binning_logN_logP():
    # plt.style.use('ggplot')
    fig = plt.figure(figsize=(12,10))
    # ax2 = fig.add_subplot(211)
    ax2 = None

    if ax2 is None:
        ax = fig.add_subplot(111)
    else:
        ax = fig.add_subplot(212, sharex=ax2)
        ax2.set_title('Comparison of different binning for logN logP from Stern+01',**font)
    # sternfile2 = 'Stern_digitized.csv'
    # mid_bins_digi = pf.read_column(sternfile2,0, splitter=', ')/0.75
    # stern_digi = pf.read_column(sternfile2,1, splitter=', ')
    #print bins
    hist_S2, bins_S2 = np.histogram(final_long_peak_flux_1024_S, bins=bins)
    hist_S2 = np.asarray(hist_S2, dtype=float)
    hist_S2_err = np.sqrt(hist_S2)

    # change the bins to avoid instabilities in Chi2
    bins2 = []
    for i in range(len(bins)-8):
        bins2.append(bins[i])
    #print bins[-8:]
    bins2.append(16.)
    bins2.append(20.)
    bins2.append(28.)
    bins2.append(50.)
    bins2 = np.asarray(bins2, dtype=float)
    bins_middle_v2 = np.zeros(len(bins2)-1)
    bins_errp_v2 = np.zeros(len(bins_middle_v2))
    bins_errm_v2 = np.zeros(len(bins_middle_v2))

    for i in range(len(bins_middle_v2)):
        bins_middle_v2[i] = np.sqrt(bins2[i]*bins2[i+1])
        bins_errp_v2[i] = bins2[i+1] - bins_middle_v2[i]
        bins_errm_v2[i] = bins_middle_v2[i] - bins2[i]

    hist_S_v2, __unused = np.histogram(final_long_peak_flux_1024_S, bins=bins2)
    ln_oi = 0.
    print(f"i, bin : N, ln(N!), sum(ln(N!)")
    for i, val in enumerate(hist_S_v2):
        ln_oi += val*np.log(val) - val
        print(f"bin {i}, {bins2[i]}: {val}, {val*np.log(val) - val}, {ln_oi}")
    print(f"Ln(o_i) = {ln_oi} before efficiency correction")
    hist_S_v2 = np.asarray(hist_S_v2, dtype=float)
    hist_S_v2_err = np.sqrt(hist_S_v2)


    if ax2 is not None:
        # Plot raw histograms
        ax2.errorbar(bins_middle,    hist_S2,   yerr=hist_S2_err, xerr=[bins_errm, bins_errp],    alpha=0.8, color='teal', fmt='.')
        ax2.errorbar(bins_middle_v2, hist_S_v2, yerr=hist_S_v2_err, xerr=[bins_errm_v2, bins_errp_v2],    alpha=0.8, color='darkorange', fmt='.')
        ax2.scatter(bins_middle,    hist_S2,    alpha=0.8, marker='o', color='teal',  s=50,     label='Stern from catalog w/ old bins')
        ax2.scatter(bins_middle_v2, hist_S_v2,  alpha=0.8, marker='s', color='darkorange',s=50, label='Stern from catalog w/ new bins')

        ax2.axhline(10, ls='--', color='k')
        ax2.text(bins_middle_v2[0], 5, r'minimum to avoid $\chi^2$ instabilities')
        ax2.set_ylim(0.5, 3000)
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.set_ylabel(r'$\Delta N$ ', **font)
        plt.setp(ax2.get_xticklabels(),visible=False)

    #ax2.hist(final_long_peak_flux_1024_S, bins=bins)

    # Correct for efficiency and useful time
    hist_S_v2 /= efficiency_correction_Stern(bins_middle_v2)
    hist_S_v2_err /= efficiency_correction_Stern(bins_middle_v2)
    ln_oi = 0.
    for i, val in enumerate(hist_S_v2):
        print(f"bin {i}: {val}")
        ln_oi += val*np.log(val) - val
    print(f"Ln(o_i) = {ln_oi} after efficiency correction")

    for i in range(len(bins_middle_v2)):
        hist_S_v2[i] /= np.log10( bins2[i+1] / bins2[i] )
        hist_S_v2[i] /= 9.1*0.67
        hist_S_v2_err[i] /= np.log10( bins2[i+1] / bins2[i] )
        hist_S_v2_err[i] /= 9.1*0.67
    
    for i in range(len(bins_middle)):
        hist_S2[i] /= np.log10( bins[i+1] / bins[i] )
        hist_S2[i] /= 9.1*0.67
        hist_S2_err[i] /= np.log10( bins[i+1] / bins[i] )
        hist_S2_err[i] /= 9.1*0.67
    hist_S2 /= efficiency_correction_Stern(bins_middle)
    hist_S2_err /= efficiency_correction_Stern(bins_middle)

    # Plot
    #ax.scatter(mid_bins_digi,  stern_digi, alpha=0.6, marker='x', color='k', s=100,   label='Stern digitized (fig. 23)')
    ax.errorbar(bins_middle,    hist_S2,   yerr=hist_S2_err, xerr=[bins_errm, bins_errp],    alpha=0.8, color='teal', fmt='.')
    ax.scatter(bins_middle,    hist_S2,     alpha=0.6, marker='o', color='teal',  s=50,     label='Stern from catalog w/ old bins')
    ax.errorbar(bins_middle_v2, hist_S_v2, yerr=hist_S_v2_err, xerr=[bins_errm_v2, bins_errp_v2],    alpha=0.8, color='darkorange', fmt='.')
    ax.scatter(bins_middle_v2, hist_S_v2,  alpha=0.6, marker='s', color='darkorange',s=50, label='Stern from catalog w/ new bins')
    #ax.errorbar(bins_middle,   hist_Stern, yerr=hist_Stern_err,color='g', fmt='.')
    #ax.scatter(bins_middle,   hist_Stern, color='g',alpha=0.6, marker='D', s=50, label='Constraint from Daigne+06')
    #plot_all_Stern(ax)
    fig.subplots_adjust(hspace=0)
    ax.set_xlabel(r'Peak flux $\mathrm{[ph\,cm^{-2}\,s^{-1}\,50-300\,keV]}$', **font)
    ax.set_ylabel(r'$\Delta N/\Delta \,(\mathrm{log\,}P) \mathrm{[GRB/yr\,in\,4\pi]}$ ', **font)#\,\mathrm{[yr^{-1}]}
    ax.set_xlim(0.05, 200)
    ax.set_ylim(0.5, 3000)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend(loc=3,scatterpoints=1, numpoints=1)

    write_to_file = False
    if write_to_file:
        outputfile_Stern_rebinned = 'Stern_lognlogp_rebinned.txt'
        f = open(outputfile_Stern_rebinned, 'w')
        f.write("# Log N Log P from Stern+01 catalog, rebinned to avoid low number (<10) bins.\n")
        f.write("# Peak count rate was converted to peak flux assuming average ratio of 0.75 (pcr/pflx = 0.75).\n")
        f.write("# The last bin ends at {}, and there are 27 bins.\n".format(bins2[-1]))
        f.write("# Nb_GRBs/year/log10(pflx) is corrected for efficiency as described in Stern+01 (Fig.22).\n")
        f.write("# 1. Peak flux (left bin edge) [ph/cm2/s]  2. log10(Nb_GRBs/year/log10(pflx))  3. 1 sigma error in log10 scale\n")
        for i in range(len(hist_S_v2)):
            hist_S_v2_err[i] = np.log10((hist_S_v2[i]+hist_S_v2_err[i])/hist_S_v2[i])
            if i == 0:
                f.write("{:e} {:e} {:e}\n".format(bins2[i],np.log10(hist_Stern[i]), np.log10(hist_Stern_err[i])))
            else:
                f.write("{:e} {:e} {:e}\n".format(bins2[i],np.log10(hist_S_v2[i]), hist_S_v2_err[i]))

        f.close()

    plt.show()#1.172680e+03 1.525687e+02
    return


def plot_Stern(run, ax, label='Model',obs_leg=1, post_process_mode=False, N_constraint=6, **kwargs):
    """
        Function that plots Stern constraint onto the provided ax.
        If post_process_mode is True, it will average with weights equal to the inverse of the Chi2,
        So that models with high Chi2 have lower weight.
    """
    file_Stern     = root_dir + run +'/Stern_constraint.dat'
    file_Stern_err = root_dir + run +'/Stern_constrainterr.dat'

    P23_Stern      = pf.read_column(file_Stern, 0)
    LogN_Stern_obs = pf.read_column(file_Stern, 2)
    P23_Stern_av          = pf.read_column(file_Stern_err, 0)
    LogN_Stern_obs_av     = pf.read_column(file_Stern_err, 3)
    LogN_Stern_obs_av_err = pf.read_column(file_Stern_err, 4)
    
    ## Small code for the errorbars to work on logscale ###
    LogN_Stern_errp = LogN_Stern_obs_av * (10**LogN_Stern_obs_av_err - 1.) 
    LogN_Stern_errm = LogN_Stern_obs_av * (1. - 10**(-LogN_Stern_obs_av_err))
    P23_Stern_err   = np.zeros(len(P23_Stern_av))
    for i in range(len(P23_Stern_av)):
        j = 2*i
        P23_Stern_err[i] = np.log10(P23_Stern[j+1]/P23_Stern_av[i])
    P23_Stern_errp = P23_Stern_av * (10**P23_Stern_err - 1.)
    P23_Stern_errm = P23_Stern_av * (1. - 10**(-P23_Stern_err))


    if obs_leg == 1 :
        ax.errorbar(P23_Stern_av, LogN_Stern_obs_av, xerr=[P23_Stern_errm,P23_Stern_errp] , yerr=[LogN_Stern_errm, LogN_Stern_errp], markersize=5, marker='s',capsize=0, label='Stern+01',fmt='.', color='k')#,color=plt.rcParams['axes.color_cycle'][0])
    else :
        ax.errorbar(P23_Stern_av, LogN_Stern_obs_av, xerr=[P23_Stern_errm,P23_Stern_errp] , yerr=[LogN_Stern_errm, LogN_Stern_errp], markersize=5, marker='s',capsize=0,fmt='.', color='k')#,color=plt.rcParams['axes.color_cycle'][0])

    if post_process_mode is False :
        LogN_Stern_mod = pf.read_column(file_Stern, 1)
        LogN_Stern_mod_av     = pf.read_column(file_Stern_err, 1)
        LogN_Stern_mod_av_err = pf.read_column(file_Stern_err, 2) 
        ax.plot(P23_Stern, LogN_Stern_mod, label=label, lw=1.5, **kwargs)
        #ax.errorbar(P23_Stern_av, LogN_Stern_mod_av, yerr=LogN_Stern_mod_av_err,fmt='.', capsize=0,color='k')#plt.rcParams['axes.color_cycle'][color])
    
    else:
        # read data
        file_Stern_post_proc = root_dir + run +'/Stern_constraint_post_proc.dat'
        data = np.genfromtxt(file_Stern_post_proc, dtype=None)
        n_lines = len(data)
        print("Post-processed {} models".format(n_lines))
        data = np.transpose(data)
        n_bins = len(data)-(N_constraint+1) #because of the size of chi2 array from f90 code
        Chi2_array = data[0]

        LogN_Stern_hist_mod_med = np.zeros(2*n_bins)
        LogN_Stern_hist_mod_std = np.zeros(2*n_bins)
        for i in range(n_bins):
            LogN_Stern_hist_mod_med[2*i] = np.average(data[i+1 + N_constraint], weights=data[0]**(-1))
            LogN_Stern_hist_mod_med[2*i+1] = LogN_Stern_hist_mod_med[2*i]
            LogN_Stern_hist_mod_std[2*i] = np.std(data[i+1 + N_constraint])
            LogN_Stern_hist_mod_std[2*i+1] = LogN_Stern_hist_mod_std[2*i]           
        # plot it
        ax.plot(P23_Stern, LogN_Stern_hist_mod_med, lw='2',label=label, **kwargs)
        ax.fill_between(P23_Stern, LogN_Stern_hist_mod_med-LogN_Stern_hist_mod_std,LogN_Stern_hist_mod_med+LogN_Stern_hist_mod_std, alpha=0.3, **kwargs)

    ax.set_title('Intensity Constraint', **font)
    ax.set_xlabel(r'Peak flux $\mathrm{[ph\,cm^{-2}\,s^{-1}\,50-300\,keV]}$', **font)
    ax.set_ylabel(r' $\Delta N/\Delta \,(\mathrm{log\,}P)\,\mathrm{[yr^{-1}]}$ ', **font)
    ax.set_xlim([0.05, 150.])
    ax.set_xscale('log')
    ax.set_yscale('log')
    leg = ax.legend(loc='best', numpoints=1)
    leg.get_frame().set_edgecolor('k')
    
    return  


def plot_all_Stern(ax):
    plot_Stern(run, ax)
    sternfile2 = '/nethome/palmerio/1ere_annee/Frederic/GRB_population_code/observational_constraints/Stern_binning/Stern_digitized.csv'
    mid_bins_digi = pf.read_column(sternfile2,0, splitter=', ')/0.75
    stern_digi = pf.read_column(sternfile2,1, splitter=', ') 


    stern22 ='/nethome/palmerio/1ere_annee/Frederic/GRB_population_code/observational_constraints/Stern_binning/Value.csv'
    mid_bins_digi2 = pf.read_column(stern22,0, splitter=', ')/0.75
    stern_digi2 = pf.read_column(stern22,1, splitter=', ') 

    stern22errp ='/nethome/palmerio/1ere_annee/Frederic/GRB_population_code/observational_constraints/Stern_binning/Errors plus.csv'
    mid_bins_digi2errp = pf.read_column(stern22errp,0, splitter=', ')/0.75
    stern_digi2errp = pf.read_column(stern22errp,1, splitter=', ') 

    stern22errm ='/nethome/palmerio/1ere_annee/Frederic/GRB_population_code/observational_constraints/Stern_binning/Errors minus.csv'
    mid_bins_digi2errm = pf.read_column(stern22errm,0, splitter=', ')/0.75
    stern_digi2errm = pf.read_column(stern22errm,1, splitter=', ') 



    ax.scatter(mid_bins_digi2errp, stern_digi2errp, marker='s', color='r')
    ax.scatter(mid_bins_digi2errm, stern_digi2errm, marker='s', color='b')
    ax.scatter(mid_bins_digi2, stern_digi2, marker='+', color='k')
    ax.scatter(mid_bins_digi, stern_digi, marker='x', color='darkorange')
    ax.errorbar(bins_middle,   hist_Stern, yerr=hist_Stern_err,color='g', fmt='.')
    ax.scatter(bins_middle,   hist_Stern, color='g',alpha=0.6, marker='D', label='Constraint from Daigne+06')

    return


decide_binning_logN_logP()
