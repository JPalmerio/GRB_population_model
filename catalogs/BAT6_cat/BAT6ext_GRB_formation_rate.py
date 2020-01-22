# To allow the importation of plotting function module anywhere
import sys
import platform
if platform.system() == 'Linux':
	sys.path.insert(0,'/nethome/palmerio/Dropbox/Plotting_GUI/Src')
elif platform.system() == 'Darwin': 
	sys.path.insert(0,'/Users/palmerio/Dropbox/Plotting_GUI/Src')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import plotting_functions as pf


root_dir = '/nethome/palmerio/1ere_annee/Frederic/catalogs/BAT6_cat/'
filename = root_dir + 'BAT6ext_GRB_formation_rate.txt'
outputfilename = root_dir + 'z_cumul_distr_Pesc16.txt'


z = pf.read_column(filename, 0) 
distr = pf.read_column(filename, 1)
distr_err = pf.read_column(filename, 2)

plt.errorbar(z, distr, yerr=distr_err, fmt='.')


plt.show()

