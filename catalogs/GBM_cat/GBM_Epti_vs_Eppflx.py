import sys
import platform

if platform.system() == 'Linux':
	sys.path.insert(0,'/nethome/palmerio/Dropbox/Plotting_GUI/Src')
elif platform.system() == 'Darwin': 
	sys.path.insert(0,'/Users/palmerio/Dropbox/Plotting_GUI/Src')

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('default')
import plotting_functions as pf
from matplotlib.transforms import blended_transform_factory
import scipy.stats as stats

"""
	This code was used to check relation between Time Integrated Peak Energy and Pflx Peak Energy
	for GBM catalog

"""

plt.style.use('presentation')
matplotlib.rc('text', usetex=True)
matplotlib.rc('font',**{'family':'serif','serif':['Palatino']})

T_min = 2.05   # secondes
P23_min = 0.01 # ph/cm2/s
filename  = 'GBM_cat_complete.txt'

# Print the content of the catalog
overhead = pf.read_overhead(filename, splitter='|', stripper='|')
print len(overhead)
for i in range(len(overhead)):
		print i, overhead[i]

name = pf.read_column(filename, 0, dtype=str, stripper='|', splitter = '|')

# T90
T90 = pf.read_column(filename, 1, stripper='|', splitter = '|')
T90_mask = np.isfinite(T90)

# pflx
peak_flux_1024 = pf.read_column(filename, 34, stripper='|', splitter = '|' )
peak_flux_1024_err = pf.read_column(filename, 35, stripper='|', splitter = '|' )
peak_flux_1024_mask = np.isfinite(peak_flux_1024) 

# Band spectrum at pflx
nd_pflx_band_ampl     = pf.read_data(filename, 79, stripper='|', splitter = '|')
nd_pflx_band_epeak    = pf.read_data(filename, 82, stripper='|', splitter = '|')
nd_pflx_band_alpha    = pf.read_data(filename, 85, stripper='|', splitter = '|')
nd_pflx_band_beta     = pf.read_data(filename, 88, stripper='|', splitter = '|')
pflx_band_redchisq    = pf.read_column(filename, 99, stripper='|', splitter = '|')
pflx_band_dof         = pf.read_column(filename, 101, stripper='|', splitter = '|')
pflx_best_fitting_model = pf.read_column(filename, 21, dtype=str, stripper='|', splitter = '|')

# Band spectrum time integrated
nd_flnc_band_ampl     = pf.read_data(filename, 9, stripper='|', splitter = '|')
nd_flnc_band_epeak    = pf.read_data(filename, 12, stripper='|', splitter = '|')
nd_flnc_band_alpha    = pf.read_data(filename, 15, stripper='|', splitter = '|')
nd_flnc_band_beta     = pf.read_data(filename, 18, stripper='|', splitter = '|')
flnc_best_fitting_model = pf.read_column(filename, 23, dtype=str, stripper='|', splitter = '|')



# Define masks
mask = T90_mask & peak_flux_1024_mask & np.isfinite(nd_flnc_band_epeak[0]) & np.isfinite(nd_pflx_band_epeak[0])
# Apply masks
nd_flnc_band_epeak_masked = pf.mask_ndarray(nd_flnc_band_epeak, mask)
nd_pflx_band_epeak_masked = pf.mask_ndarray(nd_pflx_band_epeak, mask)
peak_flux_1024_masked     = peak_flux_1024[mask]
T90_masked                = T90[mask]
print 'N GRB after finite cut : ', len(peak_flux_1024_masked)


# Define masks
mask2 = (peak_flux_1024_masked > 0.9) & (T90_masked > T_min)
# Apply masks
nd_flnc_band_epeak_masked2 = pf.mask_ndarray(nd_flnc_band_epeak_masked, mask2)
nd_pflx_band_epeak_masked2 = pf.mask_ndarray(nd_pflx_band_epeak_masked, mask2)
peak_flux_1024_masked2     = peak_flux_1024_masked[mask2]
print 'N GRB after finite cut : ', len(peak_flux_1024_masked2)
#plt.hist(np.log10(peak_flux_1024_masked2))

# sort
sorter = peak_flux_1024_masked2.argsort()
sorter = sorter[::-1]
nd_flnc_band_epeak_sorted = pf.sort_ndarray(nd_flnc_band_epeak_masked2, sorter)
nd_pflx_band_epeak_sorted = pf.sort_ndarray(nd_pflx_band_epeak_masked2, sorter)
peak_flux_1024_sorted     = peak_flux_1024_masked2[sorter]


# ratio
ratio = np.log10(nd_pflx_band_epeak_sorted[0] / nd_flnc_band_epeak_sorted[0])
fig = plt.figure(figsize=(9,7),tight_layout=True)
ax = fig.add_subplot(111)
ax.hist(ratio, bins=np.linspace(ratio.min(), ratio.max(),100), normed=True, edgecolor='k', label='Histogram of ratio')
ax.set_xlabel(r'log($\frac{E_{p,pflx}}{E_{p,ti}}$)', {'size':25})
ax.set_ylabel(r'Normalized PDF',{'size':25})
x,y = pf.asym_gaussian_pdf(mu=0.055, sigma1=0.16, sigma2=0.16, x_min=-2., x_max=2., precision=500)
ax.plot(x, y, label='Approximate gaussian')
sns.kdeplot(ratio, ax=ax, label='KDE of ratio')
ax.axvline(np.median(ratio), ls='--', c='k', label='Median ratio value: {:.3f}'.format(np.median(ratio)))
ax.legend()

q = stats.mstats.mquantiles(ratio, prob=[0.16, 0.5, 0.84])
value = q[1]
errp  = q[2]-q[1]
errm  = q[1]-q[0]
ax.text(0.05, 0.5, 'log(ratio): {:.2f} +/- ({:.2f}, {:.2f})'.format(value, errp, errm) , transform=ax.transAxes, **{'size':15})
value_lin, errp_lin, errm_lin = pf.log_to_lin(value, errp, errm)
ax.text(0.05, 0.4, 'ratio: {:.2f} +/- ({:.2f}, {:.2f})'.format(value_lin, errp_lin, errm_lin) , transform=ax.transAxes, **{'size':15})


fig = plt.figure(figsize=(8,6),tight_layout=True)
ax = fig.add_subplot(111)
#ax.set_title(r'E$_{p}$ for peak-flux vs time-integrated spectrum for GBM bursts above 0.9 ph/s/cm2')
art = pf.scatter_incomplete_ndarray(ax, nd_flnc_band_epeak_sorted, nd_pflx_band_epeak_sorted, 
									label=r'N = %d'%len(nd_flnc_band_epeak_sorted[0]),
									colormap=[np.log10(peak_flux_1024_masked2),-0.06,2.2,'viridis_r'],
									edgecolor='k', lw=0.5, alpha_errorbar=0.5, errlw=0.5)


#fit = np.polyfit(np.log10(nd_flnc_band_epeak_sorted[0]),np.log10(nd_pflx_band_epeak_sorted[0]), deg=1)
#b = fit[1]
#a = fit[0]

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'E$_{p,ti}$ [keV]')
ax.set_ylabel(r'E$_{p,pflx}$ [keV]')
x = np.logspace(-1, 7)
ax.plot(x,x,lw=2, c='k')
#y = a * x + b
#ax.plot(x,y, lw=2, c='C2', label='Linear fit with (a,b) = ({:.3f},{:.3f})'.format(a,b))

ax.legend()
ax.set_xlim(1,1e4)
ax.set_ylim(1,1e4)
cb = fig.colorbar(art)
cb.set_label('log(pflx)')



plt.show()


