# To allow the importation of plotting function module anywhere
import sys
import platform
if platform.system() == 'Linux':
	sys.path.insert(0,'/nethome/palmerio/Dropbox/Plotting_GUI/Src')
elif platform.system() == 'Darwin': 
	sys.path.insert(0,'/Users/palmerio/Dropbox/Plotting_GUI/Src')
import plotting_functions as pf

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import statsmodels.api as sm


"""
	See GBM_cat_cutoff.py for the real logNlogP figure

"""

#plt.style.use('ggplot')
homedir = '/nethome/palmerio/1ere_annee/Frederic/'
text_size = 22
font = {'fontname':'Serif', 'size':text_size}


# Swift data
filename = homedir + 'catalogs/Swift_cat/Swift_pflx_cat.txt'
Swift_name = pf.read_column(filename, 0, dtype=str)
Swift_pflx = pf.read_column(filename, 1)
Swift_pflx_err = pf.read_column(filename, 2)
Swift_pflx_masked = Swift_pflx[np.isfinite(Swift_pflx)]
Swift_pflx_err_masked = Swift_pflx_err[np.isfinite(Swift_pflx_err)]

print "min, max, median : ", min(Swift_pflx_masked), max(Swift_pflx_masked), np.median(Swift_pflx_masked)

fig = plt.figure(tight_layout=True, figsize=(10,8))
ax  = fig.add_subplot(111)
_unused, bins_logscale = np.histogram(np.log10(Swift_pflx_masked), bins=20)
ax.hist(Swift_pflx_masked, bins=10**bins_logscale, color='darkorange', label='Swift peak flux distribution')
ax.axvline(2.6, ls='--', c='k')#, label=r'$2.6\,\mathrm{ph\,cm^{-2}\,s^{-1}\,[15-150\,keV]}$')

leg=ax.legend(loc='best')
ax.set_xlabel(r'Peak flux $\mathrm{[ph\,cm^{-2}\,s^{-1}\,15-150\,keV]}$', **font)
ax.set_ylabel(r'Number of bursts', **font)
ax.set_title('Swift peak flux distribution for LGRBs', **font)
leg.get_frame().set_edgecolor('k')
ax.set_xscale('log')
#ax.annotate(r'$2.6\,\mathrm{ph\,cm^{-2}\,s^{-1}}$', xy=(2.6, 105), xycoords='data',xytext=(0, 150), textcoords='offset points',arrowprops=dict(facecolor='black', width=1, head_width=0.5),horizontalalignment='left', verticalalignment='top')

plt.show()