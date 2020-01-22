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
from scipy.stats import chi2
from scipy import interpolate
# plt.style.use('ggplot')
#plt.style.use('presentation')

matplotlib.rc('font',**{'family':'serif','serif':['Palatino']})
matplotlib.rc('text', usetex=False)

filename  = 'GBM_cat_complete.txt'
verbose = False

peak_flux_1024_nd = pf.read_data(filename, 34, stripper='|', splitter = '|' ,single_err=True)  # Batse flux (50-300 keV)
pflx_band_alpha_nd = pf.read_data(filename, 85, stripper='|', splitter = '|')
pflx_band_beta_nd  = pf.read_data(filename, 88, stripper='|', splitter = '|')

fig = plt.figure(figsize=(12,10), tight_layout=True)
axa = fig.add_subplot(211)
axb = fig.add_subplot(212, sharex=axa)
pflx_band_alpha_nd[0] = -pflx_band_alpha_nd[0]
pflx_band_beta_nd[0] = -pflx_band_beta_nd[0]
pf.scatter_incomplete_ndarray(axa,peak_flux_1024_nd,pflx_band_alpha_nd, color='k', alpha=0.9, alpha_errorbar=0.1,s=50, label='All')
pf.scatter_incomplete_ndarray(axb,peak_flux_1024_nd,pflx_band_beta_nd,  color='k', alpha=0.9, alpha_errorbar=0.1,s=50, label='All')
#axa.set_yscale('log')
axa.set_xscale('log')
axb.set_yscale('log')
#plt.show()

overhead = pf.read_overhead(filename, splitter='|', stripper='|')
print len(overhead)
for i in range(len(overhead)):
		print i, overhead[i]
#raise SystemExit


name = pf.read_column(filename, 0, dtype=str, stripper='|', splitter = '|')
T90 = pf.read_column(filename, 1, stripper='|', splitter = '|')
T90_mask = np.isfinite(T90)

peak_flux_1024 = pf.read_column(filename, 34, stripper='|', splitter = '|' )  # Batse flux (50-300 keV)
peak_flux_1024_err = pf.read_column(filename, 35, stripper='|', splitter = '|' )  # Batse flux (50-300 keV)
peak_flux_1024_mask = np.isfinite(peak_flux_1024)

T90_masked = T90[T90_mask & peak_flux_1024_mask]
peak_flux_1024_masked = peak_flux_1024[T90_mask & peak_flux_1024_mask]

print 'N_GRB avant coupure : ',len(peak_flux_1024_masked)

peak_flux_1024_timecut = peak_flux_1024_masked[T90_masked >= 2.]

print 'N_GRB apres coupure en duree : ', len(peak_flux_1024_timecut)

final_peak_flux_1024 = peak_flux_1024_timecut[peak_flux_1024_timecut >= 0.9]

print 'N_GRB apres coupure en duree et en flux : ', len(final_peak_flux_1024)
#plt.hist(np.log10(final_peak_flux_1024))

# Band
pflx_band_ampl     = pf.read_column(filename, 79, stripper='|', splitter = '|')
pflx_band_epeak    = pf.read_column(filename, 82, stripper='|', splitter = '|')
pflx_band_alpha    = pf.read_column(filename, 85, stripper='|', splitter = '|')
pflx_band_beta     = pf.read_column(filename, 88, stripper='|', splitter = '|')
pflx_band_redchisq = pf.read_column(filename, 99, stripper='|', splitter = '|')
pflx_band_dof      = pf.read_column(filename, 101, stripper='|', splitter = '|')
pflx_best_fitting_model = pf.read_column(filename, 21, dtype=str, stripper='|', splitter = '|')
# positive error
pflx_band_ampl_errorp  = pf.read_column(filename, 80, stripper='|', splitter = '|')
pflx_band_epeak_errorp = pf.read_column(filename, 83, stripper='|', splitter = '|')
pflx_band_alpha_errorp = pf.read_column(filename, 86, stripper='|', splitter = '|')
pflx_band_beta_errorp  = pf.read_column(filename, 89, stripper='|', splitter = '|')
# negative error
pflx_band_ampl_errorm  = pf.read_column(filename, 81, stripper='|', splitter = '|')
pflx_band_epeak_errorm = pf.read_column(filename, 84, stripper='|', splitter = '|')
pflx_band_alpha_errorm = pf.read_column(filename, 87, stripper='|', splitter = '|')
pflx_band_beta_errorm  = pf.read_column(filename, 90, stripper='|', splitter = '|')

# Compton
pflx_comp_ampl  = pf.read_column(filename, 58, stripper='|', splitter = '|')
pflx_comp_epeak = pf.read_column(filename, 61, stripper='|', splitter = '|')
pflx_comp_index = pf.read_column(filename, 64, stripper='|', splitter = '|')

N_plaw = 0
N_band = 0
N_comp = 0
N_sbpl = 0
for i in range(len(pflx_best_fitting_model)):
	pflx_best_fitting_model[i] = pflx_best_fitting_model[i].strip()
	if pflx_best_fitting_model[i] == 'pflx_plaw':
		N_plaw += 1
	elif pflx_best_fitting_model[i] == 'pflx_band': 
		N_band += 1
	elif pflx_best_fitting_model[i] == 'pflx_comp': 
		N_comp += 1
	elif pflx_best_fitting_model[i] == 'pflx_sbpl': 
		N_sbpl += 1



print 'Number of GRBs with plaw as best fit : ', N_plaw
print 'Number of GRBs with band as best fit : ', N_band
print 'Number of GRBs with comp as best fit : ', N_comp
print 'Number of GRBs with sbpl as best fit : ', N_sbpl
print 'Total number of GRBS : ', N_sbpl+ N_comp + N_band + N_plaw
print 'Total length of catalog : ', len(pflx_best_fitting_model)


# Masks
pflx_band_ampl_mask     = np.isfinite(pflx_band_ampl )
pflx_band_epeak_mask    = np.isfinite(pflx_band_epeak)
pflx_band_alpha_mask    = np.isfinite(pflx_band_alpha)
pflx_band_beta_mask     = np.isfinite(pflx_band_beta )
pflx_band_dof_mask      = np.isfinite(pflx_band_dof)
pflx_band_redchisq_mask = np.isfinite(pflx_band_redchisq)
# positive error masks
pflx_band_ampl_errorp_mask  = np.isfinite(pflx_band_ampl_errorp)
pflx_band_epeak_errorp_mask = np.isfinite(pflx_band_epeak_errorp)
pflx_band_alpha_errorp_mask = np.isfinite(pflx_band_alpha_errorp)
pflx_band_beta_errorp_mask  = np.isfinite(pflx_band_beta_errorp)
# negative error masks
pflx_band_ampl_errorm_mask  = np.isfinite(pflx_band_ampl_errorm)
pflx_band_epeak_errorm_mask = np.isfinite(pflx_band_epeak_errorm)
pflx_band_alpha_errorm_mask = np.isfinite(pflx_band_alpha_errorm)
pflx_band_beta_errorm_mask  = np.isfinite(pflx_band_beta_errorm)

# Compton masks
pflx_comp_ampl_mask  = np.isfinite(pflx_comp_ampl )
pflx_comp_epeak_mask = np.isfinite(pflx_comp_epeak)
pflx_comp_index_mask = np.isfinite(pflx_comp_index)


# fig2 = plt.figure()
# ax=fig2.add_subplot(111)

compare_pflx_band_alpha = pflx_band_alpha[pflx_comp_index_mask & pflx_band_alpha_mask & T90_mask]
compare_pflx_comp_index = pflx_comp_index[pflx_comp_index_mask & pflx_band_alpha_mask & T90_mask]

compare_T90 = T90[pflx_comp_index_mask & pflx_band_alpha_mask & T90_mask]

compare_pflx_BATSE = peak_flux_1024[pflx_comp_index_mask & pflx_band_alpha_mask & T90_mask]

#ax.scatter(compare_pflx_comp_index, compare_pflx_band_alpha, alpha=0.8, s=8, marker='o',  color='red', label='no cut')
#print "N_GRB with band_alpha and comp_index before cuts : ", len(compare_pflx_comp_index)

# execute time cut
compare_pflx_band_alpha_timecut = compare_pflx_band_alpha[compare_T90 >= 2]
compare_pflx_comp_index_timecut = compare_pflx_comp_index[compare_T90 >= 2]
compare_pflx_BATSE = compare_pflx_BATSE[compare_T90 >= 2.]
#print "N_GRB with band_alpha and comp_index after time cut : ", len(compare_pflx_comp_index_timecut)

# execute peak flux cut
compare_pflx_band_alpha_timecut_pflxcut = compare_pflx_band_alpha_timecut[compare_pflx_BATSE >= 0.9]
compare_pflx_comp_index_timecut_pflxcut = compare_pflx_comp_index_timecut[compare_pflx_BATSE >= 0.9]
#print "N_GRB with band_alpha and comp_index after time + pflx cut : ", len(compare_pflx_comp_index_timecut_pflxcut)


# ax.scatter(compare_pflx_comp_index_timecut_pflxcut, compare_pflx_band_alpha_timecut_pflxcut, alpha=0.8, s=3, marker='o', color='blue', label='time + pflx cut')
# x = np.linspace(-5, 6)
# ax.plot(x,x, color='k')
# ax.set_xlabel('compare_pflx_comp_index')
# ax.set_ylabel('compare_pflx_band_alpha')
# ax.legend(loc='best')


# Define sample where Band is best fit
pflx_band_ampl_best  = pflx_band_ampl[pflx_best_fitting_model == 'pflx_band']
pflx_band_epeak_best = pflx_band_epeak[pflx_best_fitting_model == 'pflx_band']
pflx_band_alpha_best = pflx_band_alpha[pflx_best_fitting_model == 'pflx_band']
pflx_band_beta_best  = pflx_band_beta[pflx_best_fitting_model == 'pflx_band']
T90_best             = T90[pflx_best_fitting_model == 'pflx_band']
pflx_BATSE_1024_best = peak_flux_1024[pflx_best_fitting_model == 'pflx_band']
# Define masks
pflx_band_ampl_best_mask  = np.isfinite(pflx_band_ampl_best)
pflx_band_epeak_best_mask = np.isfinite(pflx_band_epeak_best)
pflx_band_alpha_best_mask = np.isfinite(pflx_band_alpha_best)
pflx_band_beta_best_mask  = np.isfinite(pflx_band_beta_best)
T90_best_mask             = np.isfinite(T90_best)
pflx_BATSE_1024_best_mask = np.isfinite(pflx_BATSE_1024_best)
# Apply masks
pflx_band_ampl_best_masked  = pflx_band_ampl_best[pflx_band_ampl_best_mask]
pflx_band_epeak_best_masked = pflx_band_epeak_best[pflx_band_epeak_best_mask]
pflx_band_alpha_best_masked = pflx_band_alpha_best[pflx_band_alpha_best_mask]
pflx_band_beta_best_masked  = pflx_band_beta_best[pflx_band_beta_best_mask]
T90_best_masked             = T90_best[T90_best_mask]
pflx_BATSE_1024_best_masked = pflx_BATSE_1024_best[pflx_BATSE_1024_best_mask]
# Apply T90 cut (T90 >= 2s)
pflx_band_ampl_best_masked_timecut  = pflx_band_ampl_best_masked [T90_best_masked >= 2.]
pflx_band_epeak_best_masked_timecut = pflx_band_epeak_best_masked[T90_best_masked >= 2.]
pflx_band_alpha_best_masked_timecut = pflx_band_alpha_best_masked[T90_best_masked >= 2.]
pflx_band_beta_best_masked_timecut  = pflx_band_beta_best_masked [T90_best_masked >= 2.]
pflx_BATSE_1024_best_masked_timecut = pflx_BATSE_1024_best_masked[T90_best_masked >= 2.]
# Apply BATSE peak flux cut (pflx >= 0.9 ph/cm2/s)
pflx_band_ampl_best_masked_timecut_fluxcut  = pflx_band_ampl_best_masked_timecut [pflx_BATSE_1024_best_masked_timecut >= 0.9] #ph/cm2/s 
pflx_band_epeak_best_masked_timecut_fluxcut = pflx_band_epeak_best_masked_timecut[pflx_BATSE_1024_best_masked_timecut >= 0.9] #ph/cm2/s 
pflx_band_alpha_best_masked_timecut_fluxcut = pflx_band_alpha_best_masked_timecut[pflx_BATSE_1024_best_masked_timecut >= 0.9] #ph/cm2/s 
pflx_band_beta_best_masked_timecut_fluxcut  = pflx_band_beta_best_masked_timecut [pflx_BATSE_1024_best_masked_timecut >= 0.9] #ph/cm2/s 
pflx_BATSE_1024_best_masked_timecut_fluxcut = pflx_BATSE_1024_best_masked_timecut[pflx_BATSE_1024_best_masked_timecut >= 0.9]




#ampl
#epeak
#alpha
#beta
# Define good sample (low-energy power-law error < 0.4, high-energy power-law error < 1.0, all other parameters relative error < 0.4 from Gruber 2014)
# Apply masks
pflx_band_ampl_masked_by_errors  = pflx_band_ampl [pflx_band_ampl_errorp_mask  & pflx_band_ampl_errorm_mask  & T90_mask & peak_flux_1024_mask]
pflx_band_epeak_masked_by_errors = pflx_band_epeak[pflx_band_epeak_errorp_mask & pflx_band_epeak_errorm_mask & T90_mask & peak_flux_1024_mask]
pflx_band_alpha_masked_by_errors = pflx_band_alpha[pflx_band_alpha_errorp_mask & pflx_band_alpha_errorm_mask & T90_mask & peak_flux_1024_mask]
pflx_band_beta_masked_by_errors  = pflx_band_beta [pflx_band_beta_errorp_mask  & pflx_band_beta_errorm_mask  & T90_mask & peak_flux_1024_mask]
T90_masked_by_errors             = T90            [pflx_band_ampl_errorp_mask  & pflx_band_ampl_errorm_mask  & T90_mask & peak_flux_1024_mask]
pflx_BATSE_1024_masked_by_errors = peak_flux_1024 [pflx_band_ampl_errorp_mask  & pflx_band_ampl_errorm_mask  & T90_mask & peak_flux_1024_mask]
# positive error
pflx_band_ampl_errorp_masked_by_errors  = pflx_band_ampl_errorp [pflx_band_ampl_errorp_mask  & pflx_band_ampl_errorm_mask  & T90_mask & peak_flux_1024_mask]
pflx_band_epeak_errorp_masked_by_errors = pflx_band_epeak_errorp[pflx_band_ampl_errorp_mask  & pflx_band_ampl_errorm_mask  & T90_mask & peak_flux_1024_mask]
pflx_band_alpha_errorp_masked_by_errors = pflx_band_alpha_errorp[pflx_band_ampl_errorp_mask  & pflx_band_ampl_errorm_mask  & T90_mask & peak_flux_1024_mask]
pflx_band_beta_errorp_masked_by_errors  = pflx_band_beta_errorp [pflx_band_ampl_errorp_mask  & pflx_band_ampl_errorm_mask  & T90_mask & peak_flux_1024_mask]
# negative error
pflx_band_ampl_errorm_masked_by_errors  = pflx_band_ampl_errorm [pflx_band_ampl_errorp_mask  & pflx_band_ampl_errorm_mask  & T90_mask & peak_flux_1024_mask]
pflx_band_epeak_errorm_masked_by_errors = pflx_band_epeak_errorm[pflx_band_ampl_errorp_mask  & pflx_band_ampl_errorm_mask  & T90_mask & peak_flux_1024_mask]
pflx_band_alpha_errorm_masked_by_errors = pflx_band_alpha_errorm[pflx_band_ampl_errorp_mask  & pflx_band_ampl_errorm_mask  & T90_mask & peak_flux_1024_mask]
pflx_band_beta_errorm_masked_by_errors  = pflx_band_beta_errorm [pflx_band_ampl_errorp_mask  & pflx_band_ampl_errorm_mask  & T90_mask & peak_flux_1024_mask]
# Create lists to fill during loop
pflx_band_ampl_good  = []
pflx_band_epeak_good = []
pflx_band_alpha_good = []
pflx_band_beta_good  = []
T90_good             = []
pflx_BATSE_1024_good = []
N_a = 0
N_b = 0
N_c = 0
N_d = 0

#print "N_GRB before criteria : ", len(T90_masked_by_errors)
# Use selection criterium :
for i in range(len(T90_masked_by_errors)):
	condition_a = (pflx_band_alpha_errorp[i] <= 0.4) and (pflx_band_alpha_errorm[i] <= 0.4)
	if condition_a:
		N_a += 1
	condition_b = (pflx_band_beta_errorp[i]  <= 1.0) and (pflx_band_beta_errorm[i]  <= 1.0)
	if condition_b:
		N_b += 1
	condition_c = (pflx_band_ampl_errorp_masked_by_errors[i]/pflx_band_ampl_masked_by_errors[i]   <= 0.4) and (pflx_band_ampl_errorm_masked_by_errors[i]/pflx_band_ampl_masked_by_errors[i]   <= 0.4)
	if condition_c:
		N_c += 1
	condition_d = (pflx_band_epeak_errorp_masked_by_errors[i]/pflx_band_epeak_masked_by_errors[i] <= 0.4) and (pflx_band_epeak_errorm_masked_by_errors[i]/pflx_band_epeak_masked_by_errors[i] <= 0.4)
	if condition_d:
		N_d += 1

	if (condition_a & condition_b & condition_c & condition_d) :
		pflx_band_ampl_good.append( pflx_band_ampl_masked_by_errors[i])
		pflx_band_epeak_good.append(pflx_band_epeak_masked_by_errors[i])
		pflx_band_alpha_good.append(pflx_band_alpha_masked_by_errors[i])
		pflx_band_beta_good.append( pflx_band_beta_masked_by_errors[i])
		T90_good.append(T90_masked_by_errors[i])
		pflx_BATSE_1024_good.append(pflx_BATSE_1024_masked_by_errors[i])

if verbose:
	print "N_GRB before cuts in good sample : ", len(pflx_band_ampl_good)
	print "condition_a (on alpha) : ", N_a
	print "condition_b (on beta) : ", N_b
	print "condition_c (on ampl) : ", N_c
	print "condition_d (on epeak) : ", N_d

# Define good sample (Chi2 <=  1 sigma value for given dof)
# Apply chi2 masks
# pflx_band_ampl_masked_by_dof  = pflx_band_ampl[pflx_band_redchisq_mask & pflx_band_dof_mask]
# pflx_band_epeak_masked_by_dof = pflx_band_epeak[pflx_band_redchisq_mask & pflx_band_dof_mask]
# pflx_band_alpha_masked_by_dof = pflx_band_alpha[pflx_band_redchisq_mask & pflx_band_dof_mask]
# pflx_band_beta_masked_by_dof  = pflx_band_beta[pflx_band_redchisq_mask & pflx_band_dof_mask]
# pflx_band_redchisq_masked_by_dof = pflx_band_redchisq[pflx_band_redchisq_mask & pflx_band_dof_mask]
# pflx_band_dof_masked_by_dof   = pflx_band_dof[pflx_band_redchisq_mask & pflx_band_dof_mask]
# T90_masked_by_dof             = T90[pflx_band_redchisq_mask & pflx_band_dof_mask]
# pflx_BATSE_1024_masked_by_dof = peak_flux_1024[pflx_band_redchisq_mask & pflx_band_dof_mask]
# Create lists to fill during loop
# pflx_band_ampl_good  = []
# pflx_band_epeak_good = []
# pflx_band_alpha_good = []
# pflx_band_beta_good  = []
# T90_good             = []
# pflx_BATSE_1024_good = []
# Apply Chi2 limit
# for i in range(len(pflx_band_epeak_masked_by_dof)):
# #	chi2_limit = chi2.ppf(0.9973, pflx_band_dof_masked_by_dof[i])
# 	chi2_limit = chi2.ppf(0.6827, pflx_band_dof_masked_by_dof[i])
# 	#print "dof : ", pflx_band_dof_masked_by_dof[i], " chi2_limit : ", chi2_limit
# 	if (pflx_band_redchisq_masked_by_dof[i] * pflx_band_dof_masked_by_dof[i]) <= chi2_limit :
# 		pflx_band_ampl_good.append( pflx_band_ampl_masked_by_dof[i])
# 		pflx_band_epeak_good.append(pflx_band_epeak_masked_by_dof[i])
# 		pflx_band_alpha_good.append(pflx_band_alpha_masked_by_dof[i])
# 		pflx_band_beta_good.append( pflx_band_beta_masked_by_dof[i])
# 		T90_good.append(T90_masked_by_dof[i])
# 		pflx_BATSE_1024_good.append(pflx_BATSE_1024_masked_by_dof[i])

# Turn into array
pflx_band_ampl_good  = np.asarray(pflx_band_ampl_good )
pflx_band_epeak_good = np.asarray(pflx_band_epeak_good)
pflx_band_alpha_good = np.asarray(pflx_band_alpha_good)
pflx_band_beta_good  = np.asarray(pflx_band_beta_good )
T90_good             = np.asarray(T90_good)
pflx_BATSE_1024_good = np.asarray(pflx_BATSE_1024_good)
# Create masks for good sample
pflx_band_ampl_good_mask  = np.isfinite(pflx_band_ampl_good )
pflx_band_epeak_good_mask = np.isfinite(pflx_band_epeak_good)
pflx_band_alpha_good_mask = np.isfinite(pflx_band_alpha_good)
pflx_band_beta_good_mask  = np.isfinite(pflx_band_beta_good )
T90_good_mask             = np.isfinite(T90_good)
pflx_BATSE_1024_good_mask = np.isfinite(pflx_BATSE_1024_good)
# Apply masks to good sample
pflx_band_ampl_good_masked  = pflx_band_ampl_good [pflx_band_ampl_good_mask]
pflx_band_epeak_good_masked = pflx_band_epeak_good[pflx_band_epeak_good_mask]
pflx_band_alpha_good_masked = pflx_band_alpha_good[pflx_band_alpha_good_mask]
pflx_band_beta_good_masked  = pflx_band_beta_good [pflx_band_beta_good_mask]
T90_good_masked             = T90_good[T90_good_mask]
pflx_BATSE_1024_good_masked = pflx_BATSE_1024_good[pflx_BATSE_1024_good_mask]
# Apply T90 cut (T90 >= 2s)
pflx_band_ampl_good_masked_timecut  = pflx_band_ampl_good_masked [T90_good_masked >= 2.]
pflx_band_epeak_good_masked_timecut = pflx_band_epeak_good_masked[T90_good_masked >= 2.]
pflx_band_alpha_good_masked_timecut = pflx_band_alpha_good_masked[T90_good_masked >= 2.]
pflx_band_beta_good_masked_timecut  = pflx_band_beta_good_masked [T90_good_masked >= 2.]
pflx_BATSE_1024_good_masked_timecut = pflx_BATSE_1024_good_masked[T90_good_masked >= 2.]
# Apply BATSE peak flux cut (pflx >= 0.9 ph/cm2/s)
pflx_band_ampl_good_masked_timecut_fluxcut  = pflx_band_ampl_good_masked_timecut [pflx_BATSE_1024_good_masked_timecut >= 0.9] #ph/cm2/s 
pflx_band_epeak_good_masked_timecut_fluxcut = pflx_band_epeak_good_masked_timecut[pflx_BATSE_1024_good_masked_timecut >= 0.9] #ph/cm2/s 
pflx_band_alpha_good_masked_timecut_fluxcut = pflx_band_alpha_good_masked_timecut[pflx_BATSE_1024_good_masked_timecut >= 0.9] #ph/cm2/s 
pflx_band_beta_good_masked_timecut_fluxcut  = pflx_band_beta_good_masked_timecut [pflx_BATSE_1024_good_masked_timecut >= 0.9] #ph/cm2/s 
pflx_BATSE_1024_good_masked_timecut_fluxcut = pflx_BATSE_1024_good_masked_timecut[pflx_BATSE_1024_good_masked_timecut >= 0.9]


axa.scatter(pflx_BATSE_1024_good_masked_timecut_fluxcut, -pflx_band_alpha_good_masked_timecut_fluxcut, marker='o', color='#2ca25f', s=30, zorder=10, label='Good')
axb.scatter(pflx_BATSE_1024_good_masked_timecut_fluxcut, -pflx_band_beta_good_masked_timecut_fluxcut,  marker='o', color='#2ca25f', s=30, zorder=10, label='Good')
axa.scatter(pflx_BATSE_1024_best_masked_timecut_fluxcut, -pflx_band_alpha_best_masked_timecut_fluxcut, marker='o', color='#c51b8a', s=20, zorder=11, label='Best')
axb.scatter(pflx_BATSE_1024_best_masked_timecut_fluxcut, -pflx_band_beta_best_masked_timecut_fluxcut,  marker='o', color='#c51b8a', s=20, zorder=11, label='Best')

axa.axvline(0.9, lw=2, color='r', label='Peak flux cut-off', zorder=12)
axb.axvline(0.9, lw=2, color='r', label='Peak flux cut-off', zorder=12)
axa.set_ylim(-2.5,3)
axa.set_xlim(1e-1,1e2)
axb.set_ylim(1,50)
axa.set_ylabel(r'Band $\alpha$')
axb.set_ylabel(r'Band $\beta$')
axb.set_xlabel(r'Peak flux [$\rm ph\,s^{-1}\,cm^{-2}$]')
axa.legend()
axb.legend()
# plt.show()

# raise SystemExit

# Apply mask to total sample (regardless of best fit or chi2)
pflx_band_ampl_masked  = pflx_band_ampl [pflx_band_ampl_mask  & T90_mask & peak_flux_1024_mask]
pflx_band_epeak_masked = pflx_band_epeak[pflx_band_epeak_mask & T90_mask & peak_flux_1024_mask]
pflx_band_alpha_masked = pflx_band_alpha[pflx_band_alpha_mask & T90_mask & peak_flux_1024_mask]
pflx_band_beta_masked  = pflx_band_beta [pflx_band_beta_mask  & T90_mask & peak_flux_1024_mask]
pflx_band_redchisq_masked = pflx_band_redchisq[pflx_band_redchisq_mask& T90_mask & peak_flux_1024_mask]
T90_masked             = T90[pflx_band_ampl_mask & T90_mask & peak_flux_1024_mask]
pflx_BATSE_1024_masked = peak_flux_1024[pflx_band_ampl_mask & T90_mask & peak_flux_1024_mask]

# Apply T90 cut (T90 >= 2s)
pflx_band_ampl_masked_timecut  = pflx_band_ampl_masked [T90_masked >= 2.]
pflx_band_epeak_masked_timecut = pflx_band_epeak_masked[T90_masked >= 2.]
pflx_band_alpha_masked_timecut = pflx_band_alpha_masked[T90_masked >= 2.]
pflx_band_beta_masked_timecut  = pflx_band_beta_masked [T90_masked >= 2.]
pflx_BATSE_1024_masked_timecut = pflx_BATSE_1024_masked[T90_masked >= 2.]
# Apply BATSE peak flux cut (pflx >= 0.9 ph/cm2/s)
pflx_band_ampl_masked_timecut_fluxcut  = pflx_band_ampl_masked_timecut [pflx_BATSE_1024_masked_timecut >= 0.9] #ph/cm2/s 
pflx_band_epeak_masked_timecut_fluxcut = pflx_band_epeak_masked_timecut[pflx_BATSE_1024_masked_timecut >= 0.9] #ph/cm2/s 
pflx_band_alpha_masked_timecut_fluxcut = pflx_band_alpha_masked_timecut[pflx_BATSE_1024_masked_timecut >= 0.9] #ph/cm2/s 
pflx_band_beta_masked_timecut_fluxcut  = pflx_band_beta_masked_timecut [pflx_BATSE_1024_masked_timecut >= 0.9] #ph/cm2/s 

print 'N_GRB apres coupure en duree, en flux et finite cut : ', len(pflx_band_epeak_masked_timecut_fluxcut)

#bins_ampl  = np.linspace(0, 8*np.median(pflx_band_ampl_masked), 41)
#bins_epeak = np.logspace(np.log10(min(pflx_band_epeak_masked)), np.log10(max(pflx_band_epeak_masked)), 41)
#bins_alpha = np.linspace(min(pflx_band_alpha_masked), max(pflx_band_alpha_masked), 41)
#bins_beta  = np.linspace(min(pflx_band_beta_masked ), max(pflx_band_beta_masked ), 41)
#
#bins_ampl_best  = np.linspace(0, 8*np.median(pflx_band_ampl_best_masked), 7)
#bins_epeak_best = np.logspace(np.log10(min(pflx_band_epeak_best_masked)), np.log10(max(pflx_band_epeak_best_masked)), 12)
#bins_alpha_best = np.linspace(min(pflx_band_alpha_best_masked), max(pflx_band_alpha_best_masked), 7)
#bins_beta_best  = np.linspace(min(pflx_band_beta_best_masked ), max(pflx_band_beta_best_masked ), 7)
#	
#bins_ampl_good  = np.linspace(0, 8*np.median(pflx_band_ampl_good_masked), 21)
#bins_epeak_good = np.logspace(np.log10(min(pflx_band_epeak_good_masked)), np.log10(max(pflx_band_epeak_good_masked)), 21)
#bins_alpha_good = np.linspace(min(pflx_band_alpha_good_masked), max(pflx_band_alpha_good_masked), 41)
#bins_beta_good  = np.linspace(min(pflx_band_beta_good_masked ), max(pflx_band_beta_good_masked ), 41)


# fig = plt.figure(tight_layout=True)
# ax1 = fig.add_subplot(111)

# ax1.scatter(-pflx_band_alpha_masked_timecut_fluxcut, -pflx_band_beta_masked_timecut_fluxcut, s=40, alpha=0.6, c='gray', label='Total N=%d'%len(pflx_band_alpha_masked_timecut_fluxcut))
# ax1.scatter(-pflx_band_alpha_good_masked_timecut_fluxcut, -pflx_band_beta_good_masked_timecut_fluxcut, s=40, alpha=0.8, c='g',label='Good N=%d'%len(pflx_band_alpha_good_masked_timecut_fluxcut))
# ax1.set_xlabel(r'$\alpha_{Band}$')
# ax1.set_ylabel(r'$\beta_{Band}$')
# ax1.legend(loc='best')



bins_ampl = np.arange(0., 0.3, 0.01)
bins_epeak = np.asarray([10., 30., 60., 100., 180., 320., 565., 1000., 1800., 10000.])#np.logspace(1, 4, 13)
bins_alpha = np.arange(-1.9, 1.1, 0.2)
bins_beta = np.arange(-22., -1.9999, 0.3)
#bins_beta = np.asarray([ -12.,-8., -6., -5., -4., -3.5, -3.25, -3., -2.8, -2.6, -2.4, -2.2, -2.001])

#plt.style.use('ggplot')
fig = plt.figure(tight_layout=True)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)



ax1.hist(pflx_band_ampl_masked_timecut_fluxcut      ,bins=bins_ampl      ,normed = True, zorder=100, alpha=0.6, label=r'Total N = %d' %len(pflx_band_ampl_masked_timecut_fluxcut))
ax2.hist(pflx_band_epeak_masked_timecut_fluxcut     ,bins=bins_epeak     ,normed = True, zorder=100, alpha=0.6, label=r'Total N = %d' %len(pflx_band_epeak_masked_timecut_fluxcut))
ax3.hist(pflx_band_alpha_masked_timecut_fluxcut     ,bins=bins_alpha     ,normed = True, zorder=100, alpha=0.6, label=r'Total N = %d' %len(pflx_band_alpha_masked_timecut_fluxcut))
ax4.hist(pflx_band_beta_masked_timecut_fluxcut      ,bins=bins_beta      ,normed = True, zorder=100, alpha=0.6, label=r'Total N = %d' %len(pflx_band_beta_masked_timecut_fluxcut))
     
ax1.hist(pflx_band_ampl_best_masked_timecut_fluxcut, bins=bins_ampl  , normed = True, zorder=100, alpha=0.6, label=r'Best N = %d' %len(pflx_band_ampl_best_masked_timecut_fluxcut))
#ax2.hist(pflx_band_epeak_best_masked_timecut_fluxcut,bins=bins_epeak , zorder=100, alpha=0.6, label='Best N = %d' %len(pflx_band_epeak_best_masked_timecut_fluxcut))
ax3.hist(pflx_band_alpha_best_masked_timecut_fluxcut,bins=bins_alpha , normed = True, zorder=100, alpha=0.6, label=r'Best N = %d' %len(pflx_band_alpha_best_masked_timecut_fluxcut))
ax4.hist(pflx_band_beta_best_masked_timecut_fluxcut, bins=bins_beta  , normed = True, zorder=100, alpha=0.6, label=r'Best N = %d' %len(pflx_band_beta_best_masked_timecut_fluxcut))
     
ax1.hist(pflx_band_ampl_good_masked_timecut_fluxcut, bins=bins_ampl  , normed = True, zorder=100, alpha=0.6, label=r'Good N = %d' %len(pflx_band_ampl_good_masked_timecut_fluxcut))
ax2.hist(pflx_band_epeak_good_masked_timecut_fluxcut,bins=bins_epeak , normed = True,zorder=100, alpha=0.6, label=r'Good N = %d' %len(pflx_band_epeak_good_masked_timecut_fluxcut))
ax3.hist(pflx_band_alpha_good_masked_timecut_fluxcut,bins=bins_alpha , normed = True, zorder=100, alpha=0.6, label=r'Good N = %d' %len(pflx_band_alpha_good_masked_timecut_fluxcut))
ax4.hist(pflx_band_beta_good_masked_timecut_fluxcut, bins=bins_beta  , normed = True, zorder=100, alpha=0.6, label=r'Good N = %d' %len(pflx_band_beta_good_masked_timecut_fluxcut))
     

ax2.set_xscale('log')

ax1.set_xlabel('pflx band ampl ')
ax2.set_xlabel('pflx band epeak')
ax3.set_xlabel('pflx band alpha')
ax4.set_xlabel('pflx band beta ')

ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='best')
ax4.legend(loc='best')


#plt.style.use('default')
fig = plt.figure(tight_layout=True, figsize=(12,6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.hist(pflx_band_alpha_good_masked_timecut_fluxcut,bins=bins_alpha , color='darkorange', zorder=100, alpha=0.6, label=r'N = %d' %len(pflx_band_alpha_good_masked_timecut_fluxcut))
ax2.hist(pflx_band_beta_good_masked_timecut_fluxcut, bins=np.linspace(-5,-2,12)  , color='brown', zorder=100, alpha=0.6, label=r'N = %d' %len(pflx_band_beta_good_masked_timecut_fluxcut))
ax1.set_xlabel(r'$\alpha$', {'size':22, 'weight':'bold'})
ax2.set_xlabel(r'$\beta$', {'size':22, 'weight':'bold'})
ax1.set_ylabel(r'Number', {'size':22})
ax2.set_ylabel(r'Number', {'size':22})


ax1.legend(loc='best')
ax2.legend(loc='best')


# GBM_epeak_hist, GBM_epeak_bins = np.histogram(pflx_band_epeak_masked_timecut_fluxcut, bins = bins_epeak)
# f = open('/nethome/palmerio/1ere_annee/Frederic/GRB_population_code/observational_constraints/Ep_GBM.txt', 'w')
# f.write('# Epeak histogram from GBM catalog extracted on 27/10/16, with a cut on time ( >= 2s) and on BATSE peak flux ( >= 0.9 ph/cm2/s) \n')
# f.write('# bins | counts \n')
# f.write("# Total number of GRB : %d \n" %sum(GBM_epeak_hist))
# for i in range(len(GBM_epeak_hist)):
# 	f.write("%e \t %e \t %e \n" %(GBM_epeak_bins[i],GBM_epeak_hist[i], np.sqrt(GBM_epeak_hist[i])))
# 	f.write("%e \t %e \t %e \n" %(GBM_epeak_bins[i+1], GBM_epeak_hist[i], np.sqrt(GBM_epeak_hist[i])))
# f.close()

# bins_alpha = np.arange(-1.9, 1.025, 0.025)
# bins_beta = np.arange(-22., -1.9999, 0.1)

# GBM_alpha_hist, GBM_alpha_bins = np.histogram(pflx_band_alpha_masked_timecut_fluxcut, bins = bins_alpha)
# GBM_beta_hist, GBM_beta_bins   = np.histogram(pflx_band_beta_masked_timecut_fluxcut,  bins = bins_beta)
# FctDistr_alpha = np.zeros(len(GBM_alpha_bins))
# for i in range(len(GBM_alpha_hist)):
# 	FctDistr_alpha[i+1] = FctDistr_alpha[i] + GBM_alpha_hist[i]
# FctDistr_alpha /= sum(GBM_alpha_hist)

# FctDistr_beta = np.zeros(len(GBM_beta_bins))
# for i in range(len(GBM_beta_hist)):
# 	FctDistr_beta[i+1] = FctDistr_beta[i] + GBM_beta_hist[i]
# FctDistr_beta /= sum(GBM_beta_hist)
# #print len(GBM_beta_hist), len(GBM_beta_bins), sum(GBM_beta_hist)
# #print GBM_beta_hist, GBM_beta_bins
# #print len(GBM_alpha_hist), len(GBM_alpha_bins), sum(GBM_alpha_hist)
# #print GBM_alpha_hist, GBM_alpha_bins

# f = open('/nethome/palmerio/1ere_annee/Frederic/catalogs/GBM_cat/alpha_GBM.txt', 'w')
# f.write('# alpha histogram from GBM catalog extracted on 27/10/16 \n')
# f.write('# bins | Distribution function \n')
# f.write("# Length of this table : %d \n" %len(FctDistr_alpha))
# for i in range(len(GBM_alpha_bins)):
# 	f.write("%e \t %e \n" %(GBM_alpha_bins[i],FctDistr_alpha[i]))
# f.close()


# f = open('/nethome/palmerio/1ere_annee/Frederic/catalogs/GBM_cat/beta_GBM.txt', 'w')
# f.write('# beta histogram from GBM catalog extracted on 27/10/16 \n')
# f.write('# bins | Distribution function \n')
# f.write("# Length of this table : %d \n" %len(FctDistr_beta))
# for i in range(len(GBM_beta_bins)):
# 	f.write("%e \t %e \n" %(GBM_beta_bins[i],FctDistr_beta[i]))
# f.close()



# interp_FctDistr_alpha = interpolate.splrep(GBM_alpha_bins, FctDistr_alpha)
# interp_FctDistr_beta = interpolate.splrep(GBM_beta_bins, FctDistr_beta)

# fig2 = plt.figure()
# ax45 = fig2.add_subplot(211)
# ax46 = fig2.add_subplot(212)

# x = np.linspace(GBM_alpha_bins[0],GBM_alpha_bins[-1], 1000)
# x2 = np.linspace(GBM_beta_bins[0],GBM_beta_bins[-1], 1000)

# ax45.plot(GBM_alpha_bins, FctDistr_alpha, label="F(x)")
# ax45.plot(x, interpolate.splev(x, interp_FctDistr_alpha ), label="interpolation")

# ax46.plot(GBM_beta_bins, FctDistr_beta, label="F(x)")
# ax46.plot(x2, interpolate.splev(x2, interp_FctDistr_beta ), label="interpolation")

# ax45.legend(loc='best')
# ax46.legend(loc='best')

#bins_alpha = np.arange(-1.0, 2.0, 0.10)
#bins_beta = np.arange(2, 22.04, 0.4)
#
#GBM_alpha_hist, GBM_alpha_bins = np.histogram(-pflx_band_alpha_good_masked_timecut_fluxcut, bins = bins_alpha)
#GBM_beta_hist, GBM_beta_bins   = np.histogram(-pflx_band_beta_good_masked_timecut_fluxcut,  bins = bins_beta)
#FctDistr_alpha = np.zeros(len(GBM_alpha_bins))
#for i in range(len(GBM_alpha_hist)):
#	FctDistr_alpha[i+1] = FctDistr_alpha[i] + GBM_alpha_hist[i]
#FctDistr_alpha /= sum(GBM_alpha_hist)
#
#FctDistr_beta = np.zeros(len(GBM_beta_bins))
#for i in range(len(GBM_beta_hist)):
#	FctDistr_beta[i+1] = FctDistr_beta[i] + GBM_beta_hist[i]
#FctDistr_beta /= sum(GBM_beta_hist)
#
#f = open('/nethome/palmerio/1ere_annee/Frederic/Model_outputs/run_test_para/alpha_GBM.txt', 'w')
#f.write('# alpha histogram from GBM catalog extracted on 27/10/16 \n')
#f.write('# bins | Distribution function \n')
#f.write("# Length of this table : %d \n" %len(FctDistr_alpha))
#f.write("%e \t 0.0 \n" %(GBM_alpha_bins[0]))
#for i in range(1,len(GBM_alpha_bins)):
#	pdf = GBM_alpha_hist[i-1]/(sum(GBM_alpha_hist)*(GBM_alpha_bins[i]-GBM_alpha_bins[i-1]))
#	f.write("%e \t %e \n" %(GBM_alpha_bins[i-1], pdf))
#	f.write("%e \t %e \n" %(GBM_alpha_bins[i], pdf))
#f.close()
#
#f = open('/nethome/palmerio/1ere_annee/Frederic/Model_outputs/run_test_para/beta_GBM.txt', 'w')
#f.write('# beta histogram from GBM catalog extracted on 27/10/16 \n')
#f.write('# bins | Distribution function \n')
#f.write("# Length of this table : %d \n" %len(FctDistr_beta))
#f.write("%e \t 0.0 \n" %(GBM_beta_bins[0]))
#for i in range(1,len(GBM_beta_bins)):
#	pdf = GBM_beta_hist[i-1]/(sum(GBM_beta_hist)*(GBM_beta_bins[i]-GBM_beta_bins[i-1]))
#	f.write("%e \t %e \n" %(GBM_beta_bins[i-1], pdf))
#	f.write("%e \t %e \n" %(GBM_beta_bins[i], pdf))
#

#f.close()

#bins = np.logspace(-0.5, 1.5, 25)
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.hist(pflx_BATSE_1024_good_masked_timecut[pflx_BATSE_1024_good_masked_timecut >= 0.9], bins=bins, normed = True, zorder=100, alpha=0.6, label='Total N = %d' %len(pflx_BATSE_1024_good_masked_timecut))
#ax.hist(pflx_BATSE_1024_masked_timecut[pflx_BATSE_1024_masked_timecut >= 0.9]           ,bins=bins, normed = True, zorder=100, alpha=0.6, label='Total N = %d' %len(pflx_BATSE_1024_masked_timecut))
#ax.legend(loc='best')

plt.show()


