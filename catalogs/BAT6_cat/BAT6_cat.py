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

filename  = '/nethome/palmerio/1ere_annee/Frederic/catalogs/BAT6_cat/eBAT6_cat.txt'


overhead = pf.read_overhead(filename, splitter='|', stripper='|')
for i in range(len(overhead)):
		print i, overhead[i]


name     = pf.read_column(filename, 0, dtype=str, stripper='|', splitter = '|')
redshift = pf.read_column(filename, 1, stripper='|', splitter = '|')
Ep       = pf.read_column(filename, 13, stripper='|', splitter = '|')
Ep_err   = pf.read_column(filename, 14, stripper='|', splitter = '|')
Lum      = pf.read_column(filename, 15, stripper='|', splitter = '|')
Lum_err  = pf.read_column(filename, 16, stripper='|', splitter = '|')


redshift_mask = np.isfinite(redshift)
Ep_mask       = np.isfinite(Ep)
Lum_mask      = np.isfinite(Lum)

redshift_masked = redshift[redshift_mask]
Ep_masked       = []
for i in range(len(Ep[Ep_mask])):
	if Ep[Ep_mask][i] >=0 :
		Ep_masked.append(Ep[Ep_mask][i])
Ep_masked = np.asarray(Ep_masked)
Ep_err_masked   = np.nan_to_num(Ep_err[Ep_mask])
Lum_masked      = Lum[Lum_mask]
Lum_err_masked  = np.nan_to_num(Lum_err[Lum_mask])



# file_redshift = '/nethome/palmerio/1ere_annee/Frederic/Model_outputs/run_test_para/redshift_BAT6ext_sample.dat'

# redshift_x   = pf.read_column(file_redshift,0)
# redshift_y   = pf.read_column(file_redshift,1)

# redshift_y /= redshift_y[-1]







#plt.style.use('ggplot')
fig = plt.figure(tight_layout=True, figsize=(8,6))
ax = fig.add_subplot(111)
#ax2 = fig.add_subplot(312)
#ax3 = fig.add_subplot(313)

#bins = np.asarray([0.,0.4, 0.7, 1., 1.5, 2., 2.5, 3., 4., 6.])
ax.set_title('Redshift distribution for the BAT6 sample')
ax.hist(redshift_masked, bins=10, lw=2,color='k', histtype='step',label='N = %d'%len(redshift_masked))
#ax.plot(redshift_x, redshift_y, color='k', lw=1.5, label='MC simulation')
ax.set_xlabel('Redshift', {'size':22})
ax.set_ylabel(r'Number of GRBs', {'size':22})
ax.set_xlim([0., 6.])
ax.legend(loc='best')
#
#bins_Ep = np.logspace(1,3.5, 15)
#ax2.set_title(r'$E_{peak}$ Distribution for the BAT6 extended sample')
#ax2.hist(Ep_masked, bins=bins_Ep, color='k', label='N = %d'%len(Ep_masked))
#ax2.set_xlabel(r'Peak Energy [keV]')
#ax2.set_ylabel(r'Number')
#ax2.legend(loc='best')
#ax2.set_xscale('log')
#
#bins_Lum = np.logspace(-2,4,15)
#ax3.set_title(r'$L_{iso}$ Distribution for the BAT6 extended sample')
#ax3.hist(Lum_masked, bins=bins_Lum, color='k', label='N = %d'%len(Lum_masked))
#ax3.set_xlabel(r'Luminosity [$10^{51}\,erg/s$]')
#ax3.set_ylabel(r'Number')
#ax3.legend(loc='best')
#ax3.set_xscale('log')




plt.show()

