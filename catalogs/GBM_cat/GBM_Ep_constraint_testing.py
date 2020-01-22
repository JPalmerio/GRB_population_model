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
from matplotlib.transforms import blended_transform_factory
plt.style.use('ggplot')



fig = plt.figure()
ax = fig.add_subplot(111)
root_dir = '/nethome/palmerio/1ere_annee/Frederic/GRB_population_code/Model_outputs/'

filename = root_dir +'run_LIA/EpGBM_constraint.dat'

Ep_bins = pf.read_data(filename, 0)
Ep_hist_mod = pf.read_data(filename, 1)
Ep_hist_obs = pf.read_data(filename, 2)
x=np.linspace(1.,4., 500)
y = max(Ep_hist_obs) * pf.gaussian(x, 2.25, 0.35)
y2 = max(Ep_hist_obs) * pf.gaussian(x, 2.25, 0.375)

ep = np.linspace(1,4, 100)
ep_gauss = pf.gaussian(ep, 2.2, 0.4)*max(Ep_hist_obs)

ax.plot(Ep_bins, Ep_hist_obs, label = 'Observations')
#ax.plot(Ep_bins, Ep_hist_mod, label = 'MC simulation')
#ax.plot(ep, ep_gauss, ls=':', lw=2)
ax.plot(x,y, label='gaussian')
ax.plot(x,y2, label='gaussian2')
ax.legend(loc='best')






plt.show()


