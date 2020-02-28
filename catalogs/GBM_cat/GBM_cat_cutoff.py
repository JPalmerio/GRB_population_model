import sys
import platform
if platform.system() == 'Linux':
    sys.path.insert(0,'/nethome/palmerio/Dropbox/Plotting_GUI/Src')
elif platform.system() == 'Darwin':
    sys.path.insert(0,'/Users/palmerio/Dropbox/Plotting_GUI/Src')

# print platform.system()
# print sys.path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import plotting_functions as pf
#from astroML.plotting import hist


plt.style.use('presentation')
matplotlib.rc('font',**{'family':'serif','serif':['Palatino']})
matplotlib.rc('text', usetex=True)

"""
    This code was used to determine the cutoff peak flux for the GBM catalog.
    data taken 25/10/16.

"""
text_size = 22
font = {'fontname':'Serif', 'size':text_size}


def efficiency_correction_Stern(pflx, c_e0=0.097, nu=2.34, norm=0.7):
    """
        The efficiency function of BATSE for detecting GRBs as a function of peak flux, derived by Stern+01
        c_e0 is in [counts/s/cm2]
        pflx is in [ph/s/cm2]
    """
    c_e = pflx * 0.75  # the conversion factor from counts to pflx comes from the Stern+01 paper as well, figure 7.
    return norm * (1.0 - np.exp(-(c_e/c_e0)**2))**nu


def create_hist(filename, n_col, bins=40, ax=None, **kwargs):

    data = pf.read_column(filename, n_col, splitter='|', stripper='|')
    data_mask = np.isfinite(data)
    data = data[data_mask]
    if ax is not None:
        ax.hist(data, bins=bins, **kwargs)
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(data, bins=bins, **kwargs)
    return


T_min = 2.05   # secondes
P23_min = 0.07  # ph/cm2/s
N_bin_min = 1  # Nb of bins with flx >= 50% of pflx -> more robust than T90>2s because of uncertainties on T90 measurement (see Stern+01)
verbose = 1

#################### Stern #######################
filename_S = 'catalogs/BATSE_cat/Stern_cat.txt'
name_S = pf.read_column(filename_S, 0, dtype=str, splitter='\t')
T90_S  =  pf.read_column(filename_S, 7, splitter='\t')
N_bin_S = pf.read_column(filename_S, 8, splitter='\t')
peak_flux_1024_S = pf.read_column(filename_S, 3, splitter='\t')
peak_flux_1024_mask_S = np.isfinite(peak_flux_1024_S)  
final_peak_flux_1024_S = peak_flux_1024_S[peak_flux_1024_mask_S]
final_long_peak_flux_1024_S = peak_flux_1024_S[N_bin_S > N_bin_min]
#final_long_peak_flux_1024_S2 = peak_flux_1024_S[T90_S >= T_min]
#print 'max pflx Stern :', max(final_long_peak_flux_1024_S)
#print 'high pflx Stern :', final_long_peak_flux_1024_S[final_long_peak_flux_1024_S>30.]
if verbose >= 1 :
    print( '----- Stern catalog -----')
    print( 'Number of GRBs before N_bin cut : %d'%len(final_peak_flux_1024_S))
    print( 'Number of GRBs after  N_bin cut : %d'%len(final_long_peak_flux_1024_S))
    print( '-------------------------')
    print( '\n')

#################### GBM #######################
filename = 'catalogs/GBM_cat/GBM_cat.txt'
# Print the content of the catalog
# overhead = pf.read_overhead(filename, splitter='|', stripper='|')
# for i in range(len(overhead)):
#       print i, overhead[i]
name = pf.read_column(filename, 0, dtype=str, stripper='|', splitter = '|')
T90 = pf.read_column(filename, 1, stripper='|', splitter = '|')
T90_mask = np.isfinite(T90)
print( 'min and max name :',np.sort(name)[0],np.sort(name)[-1])

peak_flux_1024 = pf.read_column(filename, 17, stripper='|', splitter = '|' )
peak_flux_1024_err = pf.read_column(filename, 18, stripper='|', splitter = '|' )
peak_flux_1024_mask = np.isfinite(peak_flux_1024)

final_peak_flux_1024 = peak_flux_1024[peak_flux_1024_mask & T90_mask]
final_T90 = T90[peak_flux_1024_mask & T90_mask]

final_long_peak_flux_1024_nopflxcut  = final_peak_flux_1024[final_T90 >= T_min]
final_short_peak_flux_1024_nopflxcut = final_peak_flux_1024[final_T90 <  T_min]

final_long_peak_flux_1024 = final_long_peak_flux_1024_nopflxcut[final_long_peak_flux_1024_nopflxcut >= P23_min]
final_short_peak_flux_1024 = final_short_peak_flux_1024_nopflxcut[final_short_peak_flux_1024_nopflxcut >= P23_min]
if verbose >= 1 :
    print( '----- GBM catalog -----')
    print( 'Number of GRBs before T90 cut : %d'%len(final_peak_flux_1024))
    print( 'Number of GRBs after  T90 cut : %d'%len(final_long_peak_flux_1024_nopflxcut))
    print( 'Number of GRBs after  P23 cut : %d'%len(final_long_peak_flux_1024))
    print( '-------------------------')
    print( '\n')


#################### BATSE ####################
filename2 = 'catalogs/BATSE_cat/BATSE_cat.txt'
trigger_name_B = pf.read_column(filename2, 0, dtype=str, stripper='|', splitter = '|')
name_B = pf.read_column(filename2, 1, dtype=str, stripper='|', splitter = '|')
T90_B = pf.read_column(filename2, 8, stripper='|', splitter = '|')
T90_mask_B = np.isfinite(T90_B)

peak_flux_1024_B = pf.read_column(filename2, 6, stripper='|', splitter = '|' )
peak_flux_1024_err_B = pf.read_column(filename2, 7, stripper='|', splitter = '|' )
peak_flux_1024_mask_B = np.isfinite(peak_flux_1024_B)

final_peak_flux_1024_B = peak_flux_1024_B[peak_flux_1024_mask_B & T90_mask_B]
final_T90_B = T90_B[peak_flux_1024_mask_B & T90_mask_B]

final_long_peak_flux_1024_B  = final_peak_flux_1024_B[final_T90_B >= T_min]
final_short_peak_flux_1024_B = final_peak_flux_1024_B[final_T90_B  < T_min]

if verbose >= 1 :
    print( '----- BATSE catalog -----')
    print( 'Number of GRBs before T90 cut : %d'%len(final_peak_flux_1024_B))
    print( 'Number of GRBs after  T90 cut : %d'%len(final_long_peak_flux_1024_B))
    print( '-------------------------')
    print( '\n')


#################### Swift ####################
# There is already a T90 > 2s cut on this catalog
filename3 = 'catalogs/Swift_cat/Swift_pflx_cat.txt'
trigger_name_Sw= pf.read_column(filename3, 0, dtype=str)

peak_flux_1024_Sw = pf.read_column(filename3, 1)
peak_flux_1024_err_Sw = pf.read_column(filename3, 2)
peak_flux_1024_mask_Sw = np.isfinite(peak_flux_1024_Sw)

final_peak_flux_1024_Sw = peak_flux_1024_Sw[peak_flux_1024_mask_Sw]
final_name_Sw = trigger_name_Sw[peak_flux_1024_mask_Sw]

if verbose >= 1 :
    print( '----- Swift catalog -----')
    print( 'Number of GRBs : %d'%len(final_peak_flux_1024_Sw))
    print( '-------------------------')
    print( '\n')

#### Observational constraints used in fortran code
file_Stern  = 'observational_constraints/lognlogp.stern.dat'
file_Komm   = 'observational_constraints/Kommers2.dat'

# Kommers
bins_Komm = pf.read_column(file_Komm, 0, array=False)
bins_Komm.append(20.)
bins_middle_Komm = np.zeros(len(bins_Komm)-1)
bins_errp_Komm = np.zeros(len(bins_middle_Komm))
bins_errm_Komm = np.zeros(len(bins_middle_Komm))
hist_Komm = pf.read_column(file_Komm, 2)
hist_Komm_err = pf.read_column(file_Komm, 3)

# Stern
bins = pf.read_column(file_Stern, 0, array=False)
hist_Stern = 10**pf.read_column(file_Stern, 1)
hist_Stern_err = 10**pf.read_column(file_Stern, 2)
bins.append(1.8)
bins = 10**np.asarray(bins)/0.75
bins_middle = np.zeros(len(bins)-1)
bins_errp = np.zeros(len(bins_middle))
bins_errm = np.zeros(len(bins_middle))

hist_GBM, bins_GBM = np.histogram(final_long_peak_flux_1024, bins=bins)
hist_BATSE, bins_BATSE = np.histogram(final_long_peak_flux_1024_B, bins=bins)
hist_S, bins_S = np.histogram(final_long_peak_flux_1024_S, bins=bins)

hist_Sw, bins_Sw = np.histogram(final_peak_flux_1024_Sw, bins=bins)
hist_GBM = np.asarray(hist_GBM, dtype=float)
hist_BATSE = np.asarray(hist_BATSE, dtype=float)
hist_S = np.asarray(hist_S, dtype=float)
hist_Sw = np.asarray(hist_Sw, dtype=float)
N_tot_B   = sum(hist_BATSE)
N_tot_S   = sum(hist_S)
N_tot_Sw  = sum(hist_Sw)
N_tot_GBM = sum(hist_GBM)
if verbose >= 1 :
    print( 'Number of GRBs actually used in plotting : ')
    print( 'N_tot_Batse = ', int(N_tot_B  ))
    print( 'N_tot_Stern = ', int(N_tot_S  ))
    print( 'N_tot_Swift = ', int(N_tot_Sw ))
    print( 'N_tot_GBM   = ', int(N_tot_GBM))
hist_GBM_err = np.sqrt(hist_GBM)
hist_BATSE_err = np.sqrt(hist_BATSE)
hist_S_err = np.sqrt(hist_S)
hist_Sw_err = np.sqrt(hist_Sw)

for i in range(len(hist_Komm)):
    bins_middle_Komm[i] = np.sqrt(bins_Komm[i]*bins_Komm[i+1])
    bins_errp_Komm[i]   = bins_Komm[i+1] - bins_middle_Komm[i]
    bins_errm_Komm[i]   = bins_middle_Komm[i] - bins_Komm[i]
    hist_Komm[i]     /= np.log10(bins_Komm[i+1]/bins_Komm[i])
    hist_Komm_err[i] /= np.log10(bins_Komm[i+1]/bins_Komm[i])


for i in range(len(bins_middle)):
    bins_middle[i] = np.sqrt(bins[i]*bins[i+1])
    bins_errp[i] = bins[i+1] - bins_middle[i]
    bins_errm[i] = bins_middle[i] - bins[i]

    hist_GBM[i]       /= np.log10( bins[i+1] / bins[i] )
    hist_GBM_err[i]   /= np.log10( bins[i+1] / bins[i] )
    hist_BATSE[i]     /= np.log10( bins[i+1] / bins[i] )#*0.67*4*np.pi*6.54 # Omega * T_util (yr) from Antier PhD
    hist_BATSE_err[i] /= np.log10( bins[i+1] / bins[i] )#*0.67*4*np.pi*6.54
    hist_S[i]         /= np.log10( bins[i+1] / bins[i] )#*0.67*4*np.pi*6.54
    hist_S_err[i]     /= np.log10( bins[i+1] / bins[i] )#*0.67*4*np.pi*6.54
    hist_Sw[i]        /= np.log10( bins[i+1] / bins[i] )
    hist_Sw_err[i]    /= np.log10( bins[i+1] / bins[i] )

correct_efficiency = True
if correct_efficiency:
    hist_S /= efficiency_correction_Stern(bins_middle)
    hist_S_err /= efficiency_correction_Stern(bins_middle)
    # hist_Komm /= efficiency_correction_Stern(bins_middle_Komm)
    # hist_Komm_err /= efficiency_correction_Stern(bins_middle_Komm)


def chi2_parameter_adjust_logNlogP(hist, hist_err, name,
    min_pflx=1, max_pflx=30, min_param=0.0001, max_param=10,
    nb_pts=500, show_plot=True):
    """
        Computes the chi2 to find the best parameter over the range [min_pflx, max_pflx]
    """
    parametre = np.linspace(min_param, max_param, nb_pts)
    chi2 = np.zeros(len(parametre))

    for i in range(len(parametre)):
        n_bins_used = 0
        for j in range(len(hist)):
            if bins_middle[j] >= min_pflx and bins_middle[j] <= max_pflx:
                chi2[i] += (hist[j]/parametre[i] - hist_Stern[j])**2 / ((hist_err[j]/parametre[i])**2 + hist_Stern_err[j]**2)
                n_bins_used += 1

    min_chi2 = min(chi2)
    best_parametre = parametre[chi2 == min_chi2][0]

    if show_plot:
        fig34 = plt.figure(tight_layout=True)
        ax34 = fig34.add_subplot(111)
        ax34.scatter(parametre, chi2, marker='o', label='best param = %.3lf\nfor chi2 = %.1lf\nwith %d dof' %(best_parametre,min_chi2, n_bins_used-1))
        ax34.set_ylabel(r'$\chi^2$')
        ax34.set_xlabel('parameter')
        ax34.set_title(r'Best parameter adjusment for logNlogP of %s' % name)
        ax34.legend(loc='best')
    return best_parametre


adjust_catalogs = True
if adjust_catalogs:
    T_utile_GBM = 8.25  # year
    pflx_min_GBM_adjust = 0.9
    pflx_max_GBM_adjust = 50
    best_parametre_GBM = chi2_parameter_adjust_logNlogP(hist_GBM,hist_GBM_err,'GBM', min_pflx=pflx_min_GBM_adjust,max_pflx=pflx_max_GBM_adjust)
    hist_GBM /= best_parametre_GBM
    hist_GBM_err /= best_parametre_GBM
    label_GBM = r'GBM (adjusted $<\Omega>*\mathrm{T_{utile}} =$%.1e [sr yr])'%(4*np.pi*best_parametre_GBM)

    hist_S /= 9.1*0.67
    hist_S_err /= 9.1*0.67
    hist_Komm *= 4*np.pi
    hist_Komm_err *= 4*np.pi
    #hist_GBM /=  0.82*8.25*8.7/(4*np.pi)
    #hist_GBM_err /= 0.82*8.25*8.7/(4*np.pi)

    #hist_GBM /= (8.7*0.82*8.)
    #hist_GBM_err /=(8.7*0.82*8.)
    #print 4*np.pi/(8.7*0.82*8.), 1./best_parametre_GBM
    label_Swift = r'Swift [15-150keV] unadjusted'

else:
    pflx_min_GBM_adjust = 0.9
    pflx_max_GBM_adjust = 100
    label_GBM = r'GBM unadjusted'
    label_Swift = r'Swift [15-150keV] unadjusted'


def plot_logN_logP():
    # plt.style.use('ggplot')
    fig = plt.figure(tight_layout=True, figsize=(11,7))
    # ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(111)

    hist_S[0] = hist_Stern[0]

    line1, unused_1, unused_2 = ax2.errorbar(bins_middle, hist_GBM,
                                             xerr=[bins_errm, bins_errp],
                                             yerr=hist_GBM_err,
                                             marker=None, capsize=3, alpha=0.9, zorder=1, linewidth=0.6, fmt='.')
    # line2, unused_1, unused_2 = ax2.errorbar(bins_middle, hist_BATSE,
    #                                          xerr=[bins_errm, bins_errp],
    #                                          yerr=hist_BATSE_err
    #                                          marker=None, capsize=1, alpha=0.9, zorder=1, linewidth=0.6, fmt='.',
    #                                          label=r'BATSE $<\Omega>*\mathrm{T_{utile}} =$%.1e [sr yr]'%(0.67*4*np.pi*6.54))
    line3, unused_1, unused_2 = ax2.errorbar(bins_middle, hist_Stern,
                                             xerr=[bins_errm, bins_errp],
                                             yerr=hist_Stern_err,
                                             marker=None, capsize=3, alpha=0.9, zorder=1, linewidth=0.6, fmt='.')
    line4, unused_1, unused_2 = ax2.errorbar(bins_middle_Komm, hist_Komm,
                                             xerr=[bins_errm_Komm, bins_errp_Komm],
                                             yerr=hist_Komm_err,
                                             marker=None, capsize=3, alpha=0.9, zorder=1, linewidth=0.6, fmt='.')
    line5, unused_1, unused_2 = ax2.errorbar(bins_middle, hist_S,
                                             xerr=[bins_errm, bins_errp],
                                             yerr=hist_S_err,
                                             marker=None, capsize=3, alpha=0.9, zorder=1, linewidth=0.6, fmt='.')
    # line6, unused_1, unused_2 = ax2.errorbar(bins_middle, hist_Sw,
    #                                          xerr=[bins_errm, bins_errp],
    #                                          yerr=hist_Sw_err,
    #                                          marker=None, capsize=1, alpha=0.9, zorder=1, linewidth=0.6, fmt='.',
    #                                          label=label_Swift) # temps utile BATSE * ajustement

    ax2.scatter(bins_middle, hist_GBM,
                color=line1.get_color(), zorder=5, marker='s', s=40,
                label=label_GBM)
    #ax2.scatter(bins_middle, hist_BATSE,    color=line2.get_color(), zorder=5, marker='s', s=40)
    ax2.scatter(bins_middle, hist_Stern,
                color=line3.get_color(), zorder=5, marker='d', s=40,
                label='Stern+01 (value from fig. 23 in paper)')
    ax2.scatter(bins_middle_Komm, hist_Komm,
                color=line4.get_color(), zorder=5, marker='o', s=40,
                label=r'Kommers+00 (value from tab. 2 in paper multiplied by $4\pi$)')
    ax2.scatter(bins_middle, hist_S,
                color=line5.get_color(), zorder=5, marker='*', s=80,
                label=r'Stern catalog corrected for useful time (9.1 yr) and  $\Omega=0.67*4\pi$')
    #ax2.scatter(bins_middle, hist_Sw,       color=line6.get_color(), zorder=5, marker='s', s=40)
    ax2.axvline(0.9, color='k', ls='--', label=r'Peak flux cutoff for GBM, P$_{cut}$ = 0.9 ph/cm2/s')
    ax2.axvspan(pflx_min_GBM_adjust, pflx_max_GBM_adjust,
                color='gray', alpha=0.2, zorder=0, label='Range for adjustment of GBM catalog')

    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.legend(loc='lower left', numpoints=1)

    ax2.set_title('1024ms timescale in 50-300 keV band for\nLGRBs from various catalogs.', **font)
    ax2.grid()

    ax2.set_xlabel(r'Peak flux $\mathrm{[ph\,cm^{-2}\,s^{-1}\,50-300\,keV]}$', **font)
    ax2.set_ylabel(r'$\Delta N/\Delta \,(\mathrm{log\,}P) \mathrm{[GRB/yr]}$'+ r' in $4\pi$', **font) # \,\mathrm{[yr^{-1}]}
    ax2.set_xlim(0.05, 200)
    ax2.set_ylim(0.5, 3000)

    return


plot_logN_logP()



# Compare Swift with and without redshift for logNlogP
def Swift_compare():
    fig = plt.figure(tight_layout=True, figsize=(12,10))
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax)
    file_redshift ='catalogs/Swift_cat/Swift_redshift_cleaned_table.txt'
    name_redshift = pf.read_column(file_redshift, 0, dtype=str)
    redshift=pf.read_column(file_redshift, 1)
    combined_Swift_table = []
    for i in range(len(name_redshift)):
        for j in range(len(final_name_Sw)):
            if final_name_Sw[j] == name_redshift[i]:
                combined_Swift_table.append([name_redshift[i], redshift[i], final_peak_flux_1024_Sw[j]])
    
    combined_Swift_table = np.asarray(combined_Swift_table)
    combined_name = combined_Swift_table.transpose()[0]
    combined_redshift = combined_Swift_table.transpose()[1]
    combined_pflx = np.asarray(combined_Swift_table.transpose()[2],dtype=float)
    
    hist_Sw_comb, bins_Sw_comb = np.histogram(combined_pflx, bins=bins)
    hist_Sw_comb = np.asarray(hist_Sw_comb, dtype=float)
    hist_Sw_comb_err = np.sqrt(hist_Sw_comb)
    N_tot_Sw_comb = sum(hist_Sw_comb)
    
    #best_parametre_Sw = chi2_parameter_adjust_logNlogP(hist_Sw,hist_Sw_err,'Swift',min_param=-10, max_param=100, show_plot=True)
    #label_Swift = r"Swift [15-150keV] (adjusted : $<\Omega>*\mathrm{T_{utile}} =$%.1e [sr yr])" %(4*np.pi*best_parametre_Sw)

    for i in range(len(bins_middle)):
        hist_Sw_comb[i]     /= np.log10( bins[i+1] / bins[i] )
        hist_Sw_comb_err[i] /= np.log10( bins[i+1] / bins[i] )
    #hist_Sw_comb     /= best_parametre_Sw
    #hist_Sw_comb_err /= best_parametre_Sw
    
    
    ax.errorbar(bins_middle, hist_Sw,      xerr=[bins_errm, bins_errp], yerr=hist_Sw_err,      marker=None, fmt='.', label='Swift N=%d'%int(N_tot_Sw))
    ax.errorbar(bins_middle, hist_Sw_comb, xerr=[bins_errm, bins_errp], yerr=hist_Sw_comb_err, marker=None, fmt='.', label='Swift with redshift N=%d'%int(N_tot_Sw_comb))
    
    
    ax.set_yscale('log')
    ax2.set_xscale('log')
    ax.legend(loc='best', numpoints=1)
    ax.set_title('1024ms timescale in 15-150 keV band for LGRB Swift catalogs\nuncorrected for useful time.', **font)
    ax2.set_xlabel(r'Peak flux $\mathrm{[ph\,cm^{-2}\,s^{-1}\,15-150\,keV]}$', **font)
    ax.set_ylabel(r'$\Delta N/\Delta \,(\mathrm{log\,}P)$ ', **font)
    ax2.set_ylabel('Fraction of GRBs\ndetected with redshift', **font)
    ax.grid()
    ax2.grid()

    redshift_efficiency = hist_Sw_comb/hist_Sw
    redshift_efficiency[np.where(~np.isfinite(redshift_efficiency))]=0
    redshift_efficiency_errp = np.zeros(len(redshift_efficiency))
    redshift_efficiency_errm = np.zeros(len(redshift_efficiency))
    for i in range(len(redshift_efficiency)):
        if hist_Sw_comb[i] > 0 and hist_Sw[i] > 0:
            redshift_efficiency_errp[i] = np.sqrt((hist_Sw_comb_err[i]/hist_Sw_comb[i])**2 + (hist_Sw_err[i]/hist_Sw[i])**2 )
            redshift_efficiency_errm[i] = redshift_efficiency_errp[i]
        else :
            redshift_efficiency_errp[i] = 0.0
            redshift_efficiency_errm[i] = 0.0
        if redshift_efficiency[i] + redshift_efficiency_errp[i] >= 1.0:
            redshift_efficiency_errp[i] = 1.0 - redshift_efficiency[i]
        if redshift_efficiency[i] - redshift_efficiency_errm[i] <= 0.0:
            redshift_efficiency_errm[i] += (redshift_efficiency[i] - redshift_efficiency_errm[i])
    
    ax2.fill_between(bins_middle, redshift_efficiency+redshift_efficiency_errp,redshift_efficiency-redshift_efficiency_errm, color='lightgray' )
    ax2.errorbar(bins_middle, redshift_efficiency, yerr=[redshift_efficiency_errm,redshift_efficiency_errp], marker=None, fmt='.', color='k')
    ax2.plot(bins_middle, redshift_efficiency, color='k', label='Redshift detection efficiency')
    ax2.axhline(N_tot_Sw_comb/N_tot_Sw, color='r', ls='--', lw=1.2, label='Average detected fraction : %.2lf'%(N_tot_Sw_comb/N_tot_Sw))
    ax2.legend(loc='best', numpoints=1)
    ax2.set_ylim(-0.1, 1.5)
    
    plt.setp(ax.get_xticklabels(),visible=False)
    fig.subplots_adjust(hspace=0,wspace=0)
    return

##Swift_compare()


plt.show()
