# To allow the importation of plotting function module anywhere
import sys
import platform
if platform.system() == 'Linux':
	sys.path.insert(0,'/nethome/palmerio/Dropbox/Plotting_GUI/Src')
elif platform.system() == 'Darwin': 
	sys.path.insert(0,'/Users/palmerio/Dropbox/Plotting_GUI/Src')
import plotting_functions as pf

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy as sp
#from astroML.density_estimation import EmpiricalDistribution

#from astroML.plotting import hist
#import seaborn as sns

"""
	This code is to provide info on the extended BAT6 sample : eBAT6
	The data comes from Pescalli et al. 2016, stored in the file eBAT6_cat.txt
	Note : Luminosity and peak energy are both time integrated over the duration of the burst.
	Such as : 
		1) Redshift distribution (differential and cumulative)
		2) Peak energy distribution (differential and cumulative)
		3) Luminosity distribution (differential and cumulative)
		4) Luminosity versus redshift
		5) Peak energy versus redshift
		6) Luminosity versus peak energy
"""

# Parameters of the script
verbosity = 0 # 0, 1 or 2
show_plot = True
text_size = 22
font = {'fontname':'Serif', 'size':text_size}
plt.style.use('ggplot')
np.random.seed(10)
root_dir = '/nethome/palmerio/1ere_annee/Frederic/catalogs/'
file_eBAT6 = root_dir + 'BAT6_cat/eBAT6_cat.txt'



# Beginning of script : 
if verbosity == 2:
	header = pf.read_overhead(file_eBAT6, stripper='|', splitter = '|')
	print "Header of %s : " %file_eBAT6
	for i in range(len(header)):
		print i, header[i]


### Data ###
first_col  = 1
name       = pf.read_data(file_eBAT6, first_col-1, dtype=str, err=False, stripper='|', splitter = '|')
redshift   = pf.read_data(file_eBAT6, first_col, err=False, stripper='|', splitter = '|')
alpha      = pf.read_data(file_eBAT6, first_col+1,  stripper='|', splitter = '|')
beta       = pf.read_data(file_eBAT6, first_col+4,  stripper='|', splitter = '|')
pflx       = pf.read_data(file_eBAT6, first_col+9,  stripper='|', splitter = '|', single_err=True)
Ep         = pf.read_data(file_eBAT6, first_col+12, stripper='|', splitter = '|', single_err=True)
Luminosity = pf.read_data(file_eBAT6, first_col+14, stripper='|', splitter = '|', single_err=True)


# 1) Redshift distribution
fig = plt.figure(tight_layout=True, figsize=(11,12))
ax  = fig.add_subplot(211)
ax2 = fig.add_subplot(212, sharex=ax)
#fig.subplots_adjust(hspace=0,wspace=0)

redshift_masked = redshift[0][np.isfinite(redshift[0])]
redshift_ECDF, ECDF = pf.unbinned_empirical_cdf(redshift_masked)
ax.plot(redshift_ECDF, ECDF, drawstyle='steps-post', lw=2, color='purple', label='eBAT6 Empirical cumulative distribution function\n < z > = %.2lf' %np.median(redshift_masked))
# test sampling from empirical distribution
#test_z = EmpiricalDistribution(redshift_masked).rvs(1000)
#redshift_ECDF_test, ECDF_test = pf.unbinned_empirical_cdf(test_z)
#ax.plot(redshift_ECDF_test, ECDF_test, drawstyle='steps-post', lw=2, ls='--', color='purple', label='eBAT6 Empirical cumulative distribution function\n < z > = %.2lf' %np.median(test_z))
ax2.hist(redshift_masked, bins=12, color='purple', edgecolor='k', alpha=0.8, label='eBAT6 GRB distribution\n N = %d, %.0lf%% complete' %(len(redshift_masked),float(100*len(redshift_masked)/len(redshift[0]))))
plt.setp(ax.get_xticklabels(),visible=False)
ax.set_title(r'eBAT6 redshift distribution', **font)
ax.set_ylabel(r'ECDF', **font)
ax2.set_ylabel(r'Number of GRBs', **font)
ax2.set_xlabel(r'Redshift (z)', **font)
ax.legend(loc='best')
ax2.legend(loc='best')


# 2) Peak energy distribution
fig = plt.figure(tight_layout=True, figsize=(11,12))
ax  = fig.add_subplot(211)
ax2 = fig.add_subplot(212, sharex=ax)
#fig.subplots_adjust(hspace=0,wspace=0)

# mask Ep from negative values (without redshift)
Ep_masked = []
Ep_masked2 = []

for Ep_i in Ep[0]:
	if np.isfinite(Ep_i) :
		if (Ep_i > 0) :
			Ep_masked.append(Ep_i)
			Ep_masked2.append(Ep_i)
		else :
			Ep_masked2.append(Ep_i)
Ep_masked2 = np.asarray(Ep_masked2)
neg_Ep = np.where(Ep_masked2 <= 0)                                  # find where there are negative Eps (i.e. without redshift)
#fake_z = EmpiricalDistribution(redshift_masked).rvs(len(neg_Ep[0])) # draw from empirical redshift distribution
fake_z = 2.2
Ep_masked2[neg_Ep] = -Ep_masked2[neg_Ep]/fake_z                     # estimate Ep intrinsic
Ep_masked2 = np.log10(Ep_masked2)
Ep_masked = np.log10(np.asarray(Ep_masked))

Ep_ECDF, ECDF = pf.unbinned_empirical_cdf(Ep_masked)
Ep_ECDF2, ECDF2 = pf.unbinned_empirical_cdf(Ep_masked2)
ax.plot(Ep_ECDF2, ECDF2, drawstyle='steps-post', lw=2, color='red'       , label='eBAT6 Empirical cumulative distribution function, simulating missing data.\n < Ep > = %.2lf keV' %10**np.median(Ep_masked2))
ax.plot(Ep_ECDF , ECDF , drawstyle='steps-post', lw=2, color='darkorange', label='eBAT6 Empirical cumulative distribution function\n < Ep > = %.2lf keV' %10**np.median(Ep_masked))
bins = ax2.hist(Ep_masked , bins=12, color='darkorange', edgecolor='k', alpha=1, zorder=10, label='eBAT6 GRB distribution\n N = %d, %.0lf%% complete' %(len(Ep_masked),float(100*len(Ep_masked)/len(Ep[0]))))
ax2.hist(Ep_masked2, bins=bins[1], color='red'       , edgecolor='k', alpha=1, zorder=5, label='eBAT6 GRB distribution, simulating missing data.\n N = %d, %.0lf%% complete' %(len(Ep_masked2),float(100*len(Ep_masked2)/len(Ep[0]))))
plt.setp(ax.get_xticklabels(),visible=False)
ax.set_title(r'eBAT6 peak energy distribution', **font)
ax.set_ylabel(r'ECDF', **font)
ax2.set_ylabel(r'Number of GRBs', **font)
ax2.set_xlabel(r'log ( Peak energy $\mathrm{keV} )$', **font)
ax.legend(loc='best')
ax2.legend(loc='best')


# 3) Luminosity distribution
fig = plt.figure(tight_layout=True, figsize=(11,12))
ax  = fig.add_subplot(211)
ax2 = fig.add_subplot(212, sharex=ax)
#fig.subplots_adjust(hspace=0,wspace=0)

log_Luminosity_masked = np.log10(Luminosity[0][np.isfinite(Luminosity[0])])
log_Luminosity_ECDF, ECDF = pf.unbinned_empirical_cdf(log_Luminosity_masked)
ax.plot(log_Luminosity_ECDF, ECDF, drawstyle='steps-post', lw=2, color='teal', label='eBAT6 Empirical cumulative distribution function\n < log(L) > = %.2lf' %np.median(log_Luminosity_masked))
ax2.hist(log_Luminosity_masked, bins=12, color='teal', edgecolor='k', alpha=0.8, label='eBAT6 GRB distribution\n N = %d, %.0lf%% complete' %(len(log_Luminosity_masked),float(100*len(log_Luminosity_masked)/len(Luminosity[0]))))
plt.setp(ax.get_xticklabels(),visible=False)
ax.set_title(r'eBAT6 Luminosity distribution', **font)
ax.set_ylabel(r'ECDF', **font)
ax2.set_ylabel(r'Number of GRBs', **font)
ax2.set_xlabel(r'$\mathrm{log\left(\,\frac{Luminosity}{10^{51}\,erg/s}\right)\,}$', **font)
ax.legend(loc='best')
ax2.legend(loc='best')



# 4) Luminosity versus redshift and 5) Peak energy versus redshift
fig = plt.figure(tight_layout=True, figsize=(11,14))
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212, sharex=ax)
pf.scatter_incomplete_ndarray(ax, redshift, Luminosity, y_is_log=True, c='g', marker='o', s=80, edgecolor='k', label='Luminosity for eBAT6 sample')
ax.set_title(r'Luminosity vs redshift for the eBAT6 sample', **font)
ax.set_ylabel(r'Luminosity $\mathrm{10^{51}\,erg/s}$', **font)
ax.set_yscale('log')
ax.legend(loc='best', scatterpoints=1)
plt.setp(ax.get_xticklabels(),visible=False)
pf.scatter_incomplete_ndarray(ax2, redshift, Ep, y_is_log=True, c='mediumvioletred', marker='o', s=80, edgecolor='k', label='Peak energy for eBAT6 sample')
ax2.set_title(r'Peak energy vs redshift for the eBAT6 sample', **font)
ax2.set_ylabel(r'Peak energy $\mathrm{keV}$', **font)
ax2.set_yscale('log')
ax2.set_xlabel(r'Redshift (z)', **font)
ax2.legend(loc='best', scatterpoints=1)
ax.set_xlim(0,7)


# 6) Luminosity versus peak energy
homedir = '/nethome/palmerio/1ere_annee/Frederic/'
filename=homedir+'GRB_population_code/observational_constraints/Yonetoku_04_L_vs_Ep.csv'
Ep_Y = pf.read_column(filename,0, splitter='\t|')
L_Y = pf.read_column(filename, 1, splitter='\t|')

fig = plt.figure(tight_layout=True, figsize=(10,8))
ax = fig.add_subplot(111)
ax.plot(Ep_Y, L_Y, label='Yonetoku+04', color='k', lw=2, zorder=0)
pf.scatter_incomplete_ndarray(ax, Ep, Luminosity, y_is_log=True, x_is_log=True, c='crimson', marker='o', s=80, edgecolor='k', label='eBAT6 sample')
ax.set_title(r'Luminosity vs peak energy for the eBAT6 sample', **font)
ax.set_xlabel(r'Peak energy $\mathrm{keV}$', **font)
ax.set_ylabel(r'Luminosity $\mathrm{10^{51}\,erg/s}$', **font)
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend(loc='best', scatterpoints=1)

if show_plot:
	plt.show()
