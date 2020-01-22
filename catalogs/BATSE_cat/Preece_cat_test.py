import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import platform
if platform.system() == 'Linux':
	sys.path.insert(0,'/nethome/palmerio/Dropbox/Plotting_GUI/Src')
elif platform.system() == 'Darwin': 
	sys.path.insert(0,'/Users/palmerio/Dropbox/Plotting_GUI/Src')
import plotting_functions as pf
from scipy import stats
from matplotlib.transforms import blended_transform_factory
import seaborn as sns
import pandas as pd



plt.style.use('default')
plt.style.use('ggplot')

root_dir = '/nethome/palmerio/1ere_annee/Frederic/GRB_population_code/observational_constraints/'
text_size = 22
font = {'fontname':'Serif', 'size':text_size}
colors=pf.generate_colors()

file1 = root_dir+'Preece_histo.ep.dat'
file2 = root_dir+'Preece_histo.eb.dat'
file_og = root_dir + 'preece.eb.dat'

Ep_og = pf.read_column(file_og, 0)
pk_og = pf.read_column(file_og, 1)

Ep = pf.read_column(file1, 0)
pk1_Ep = pf.read_column(file1, 1)
pk2_Ep = pf.read_column(file1, 2)
pk3_Ep = pf.read_column(file1, 3)
pk4_Ep = pf.read_column(file1, 4)

Eb = pf.read_column(file2, 0)
pk1_Eb = pf.read_column(file2, 1)
pk2_Eb = pf.read_column(file2, 2)
pk3_Eb = pf.read_column(file2, 3)
pk4_Eb = pf.read_column(file2, 4)

fig = plt.figure(tight_layout=True, figsize=(10,10))
ax  = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

#ax.plot(Ep, pk1_Ep,color=colors[0],lw=2, label='Ep pk1 all')
#ax.plot(Ep, pk2_Ep,color=colors[1],lw=2, label='Ep pk1 Band')
#ax2.plot(Eb, pk1_Eb,color=colors[0],lw=2, label='Eb pk1 all')
#ax2.plot(Eb, pk2_Eb,color=colors[1],lw=2, label='Eb pk1 Band')
# 
ax.plot(Ep, pk3_Ep,color=colors[2],lw=2, label='Ep pk2 all')
ax.plot(Ep, pk4_Ep,color=colors[3],lw=2, label='Ep pk2 Band')
ax.plot(Ep_og, pk_og,color=colors[0],lw=2, label='Hist used in model')
#
ax2.plot(Eb, pk3_Eb,color=colors[2],lw=2, label='Eb pk2 all')
ax2.plot(Eb, pk4_Eb,color=colors[3],lw=2, label='Eb pk2 Band')
ax2.plot(Ep_og, pk_og,color=colors[0],lw=2, label='Hist used in model')


ax.set_xlabel(r'Peak energy keV', **font)
ax.set_xscale('log')
ax.legend(loc='best')

ax2.set_xlabel(r'Break energy keV', **font)
ax2.set_xscale('log')
ax2.legend(loc='best')


fig = plt.figure(tight_layout=True, figsize=(10,10))
ax  = fig.add_subplot(111)
ax.plot(Ep, pk2_Ep,color=colors[1],lw=2, label='Ep pk1 Band')
ax.plot(Ep, pk4_Ep,color=colors[3],lw=2, label='Ep pk2 Band')
ax.set_xlabel(r'Peak energy keV', **font)
ax.set_xscale('log')
ax.legend(loc='best')


plt.show()

