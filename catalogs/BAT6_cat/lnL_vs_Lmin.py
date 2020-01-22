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
import pandas as pd
import scipy.integrate as integrate
import scipy.stats as stats



def extract_quantity(array_from_binary, quantity_index, model_index, subindex=None):
	"""
		Extracts a quantity from a binary-extracted array.
		Binary file is an output of the fortran MC GRB population code.
		
		Parameters:
		-----------

		quantity_index [int]
			The index where to find the property
		
		model_index [int]
			The index of the model you want to extract.
	"""
	if subindex is None:
		quantity = array_from_binary[model_index][quantity_index]
	else:
		quantity = array_from_binary[model_index][quantity_index][subindex]


	return quantity

font = {'size':20}

store_dir = '/home/versailles1NS/palmerio/'
root_dir = '/nethome/palmerio/1ere_annee/Frederic/GRB_population_code/Model_outputs/'
#root_dir = store_dir
outputdir = '/nethome/palmerio/1ere_annee/Frederic/plots/lnL_vs_params/'


param = 'Lbreak'

def plot_lnL_vs_Lmin(ax, run):
	binary_file = root_dir + run+'/reprise.dat'
	f = open(binary_file,'rb')
	dt = np.dtype("i4, 24i4, 10f8, 10f8, 6f8, 6f8, 7f8, 7f8, 5f8, 2i4, f8, 10f8, 10f8, 6f8, 500f8, 400f8, 200f8, 200f8, 360f8, 25f8, 10f8, 27f8, 9f8, 15f8, i4")
	sim_output_data = np.fromfile(f, dtype=dt, count=-1)
	f.close()


	lnL = np.zeros((len(sim_output_data),7))
	Param_Lum = np.zeros((len(sim_output_data), 10))

	for i in range(len(sim_output_data)):
		lnL[i] = extract_quantity(sim_output_data, 7, i)
		Param_Lum[i] = extract_quantity(sim_output_data, 2, i)


	
	ax.plot(np.log10(Param_Lum[:,1]), lnL[:,0]-np.max(lnL[:,0]), label='Total', marker='o', markersize=3)
	ax.plot(np.log10(Param_Lum[:,1]), lnL[:,3]-np.max(lnL[:,3]), label='Stern', marker='o', markersize=3)
	ax.plot(np.log10(Param_Lum[:,1]), lnL[:,5]-np.max(lnL[:,5]), label='EpGBM', marker='o', markersize=3)
	ax.plot(np.log10(Param_Lum[:,1]), 10*(lnL[:,6]-np.max(lnL[:,6])), label='eBAT6', marker='o', markersize=3)
	ax.set_ylabel(r'$\mathrm{ln(\mathcal{L}) - ln(\mathcal{L}_{max})}$', **font)
	ax.set_xlabel(r'$\mathrm{log(L_{min})}$',**font)
	ax.legend()
	
	return

def plot_lnL_vs_Lmax(ax, run):
	binary_file = root_dir + run+'/reprise.dat'
	f = open(binary_file,'rb')
	dt = np.dtype("i4, 24i4, 10f8, 10f8, 6f8, 6f8, 7f8, 7f8, 5f8, 2i4, f8, 10f8, 10f8, 6f8, 500f8, 400f8, 200f8, 200f8, 360f8, 25f8, 10f8, 27f8, 9f8, 15f8, i4")
	sim_output_data = np.fromfile(f, dtype=dt, count=-1)
	f.close()


	lnL = np.zeros((len(sim_output_data),7))
	Param_Lum = np.zeros((len(sim_output_data), 10))

	for i in range(len(sim_output_data)):
		lnL[i] = extract_quantity(sim_output_data, 7, i)
		Param_Lum[i] = extract_quantity(sim_output_data, 2, i)


	ax.plot(np.log10(Param_Lum[:,2]), lnL[:,0]-np.max(lnL[:,0]), label='Total', marker='o', markersize=3)
	ax.plot(np.log10(Param_Lum[:,2]), lnL[:,3]-np.max(lnL[:,3]), label='Stern', marker='o', markersize=3)
	ax.plot(np.log10(Param_Lum[:,2]), lnL[:,5]-np.max(lnL[:,5]), label='EpGBM', marker='o', markersize=3)
	ax.plot(np.log10(Param_Lum[:,2]), 10*(lnL[:,6]-np.max(lnL[:,6])), label='eBAT6', marker='o', markersize=3)
	ax.set_ylabel(r'$\mathrm{ln(\mathcal{L}) - ln(\mathcal{L}_{max})}$', **font)
	ax.set_xlabel(r'$\mathrm{log(L_{max})}$',**font)
	ax.legend()

	return


def plot_lnL_vs_Lbreak(ax, run):
	binary_file = root_dir + run+'/reprise.dat'
	f = open(binary_file,'rb')
	dt = np.dtype("i4, 24i4, 10f8, 10f8, 6f8, 6f8, 7f8, 7f8, 5f8, 2i4, f8, 10f8, 10f8, 6f8, 500f8, 400f8, 200f8, 200f8, 360f8, 25f8, 10f8, 27f8, 9f8, 15f8, i4")
	sim_output_data = np.fromfile(f, dtype=dt, count=-1)
	f.close()


	lnL = np.zeros((len(sim_output_data),7))
	Param_Lum = np.zeros((len(sim_output_data), 10))

	for i in range(len(sim_output_data)):
		lnL[i] = extract_quantity(sim_output_data, 7, i)
		Param_Lum[i] = extract_quantity(sim_output_data, 2, i)


	ax.plot(np.log10(Param_Lum[:,3]), lnL[:,0]-np.max(lnL[:,0]), label='Total', marker='o', markersize=3)
	ax.plot(np.log10(Param_Lum[:,3]), lnL[:,3]-np.max(lnL[:,3]), label='Stern', marker='o', markersize=3)
	ax.plot(np.log10(Param_Lum[:,3]), lnL[:,5]-np.max(lnL[:,5]), label='EpGBM', marker='o', markersize=3)
	ax.plot(np.log10(Param_Lum[:,3]), 10*(lnL[:,6]-np.max(lnL[:,6])), label='eBAT6', marker='o', markersize=3)
	ax.set_ylabel(r'$\mathrm{ln(\mathcal{L}) - ln(\mathcal{L}_{max})}$', **font)
	ax.set_xlabel(r'$\mathrm{log(L_{break})}$',**font)
	ax.legend()

	return

def plot_lnL_vs_slope(ax, run):
	binary_file = root_dir + run+'/reprise.dat'
	f = open(binary_file,'rb')
	dt = np.dtype("i4, 24i4, 10f8, 10f8, 6f8, 6f8, 7f8, 7f8, 5f8, 2i4, f8, 10f8, 10f8, 6f8, 500f8, 400f8, 200f8, 200f8, 360f8, 25f8, 10f8, 27f8, 9f8, 15f8, i4")
	sim_output_data = np.fromfile(f, dtype=dt, count=-1)
	f.close()


	lnL = np.zeros((len(sim_output_data),7))
	Param_Lum = np.zeros((len(sim_output_data), 10))

	for i in range(len(sim_output_data)):
		lnL[i] = extract_quantity(sim_output_data, 7, i)
		Param_Lum[i] = extract_quantity(sim_output_data, 2, i)


	ax.plot(Param_Lum[:,4], lnL[:,0]-np.max(lnL[:,0]), label='Total', marker='o', markersize=3)
	ax.plot(Param_Lum[:,4], lnL[:,3]-np.max(lnL[:,3]), label='Stern', marker='o', markersize=3)
	ax.plot(Param_Lum[:,4], lnL[:,5]-np.max(lnL[:,5]), label='EpGBM', marker='o', markersize=3)
	ax.plot(Param_Lum[:,4], 10*(lnL[:,6]-np.max(lnL[:,6])), label='eBAT6', marker='o', markersize=3)
	ax.set_ylabel(r'$\mathrm{ln(\mathcal{L}) - ln(\mathcal{L}_{max})}$', **font)
	ax.set_xlabel('Slope', **font)
	ax.legend()

	return

def plot_lnL_vs_Ep0(ax, run):
	binary_file = root_dir + run+'/reprise.dat'
	f = open(binary_file,'rb')
	dt = np.dtype("i4, 24i4, 10f8, 10f8, 6f8, 6f8, 7f8, 7f8, 5f8, 2i4, f8, 10f8, 10f8, 6f8, 500f8, 400f8, 200f8, 200f8, 360f8, 25f8, 10f8, 27f8, 9f8, 15f8, i4")
	sim_output_data = np.fromfile(f, dtype=dt, count=-1)
	f.close()


	lnL = np.zeros((len(sim_output_data),7))
	Param_Ep = np.zeros((len(sim_output_data), 6))

	for i in range(len(sim_output_data)):
		lnL[i] = extract_quantity(sim_output_data, 7, i)
		Param_Ep[i] = extract_quantity(sim_output_data, 5, i)


	ax.plot(Param_Ep[:,0], lnL[:,0]-np.max(lnL[:,0]), label='Total', marker='o', markersize=3)
	ax.plot(Param_Ep[:,0], lnL[:,3]-np.max(lnL[:,3]), label='Stern', marker='o', markersize=3)
	ax.plot(Param_Ep[:,0], lnL[:,5]-np.max(lnL[:,5]), label='EpGBM', marker='o', markersize=3)
	ax.plot(Param_Ep[:,0], 10*(lnL[:,6]-np.max(lnL[:,6])), label='eBAT6', marker='o', markersize=3)
	ax.set_ylabel(r'$\mathrm{ln(\mathcal{L}) - ln(\mathcal{L}_{max})}$', **font)
	ax.set_xlabel(r'$\mathrm{Ep_0}$', **font)
	ax.set_xscale('log')
	ax.legend()

	return

def plot_lnL_vs_sigmaL(ax, run):
	binary_file = root_dir + run+'/reprise.dat'
	f = open(binary_file,'rb')
	dt = np.dtype("i4, 24i4, 10f8, 10f8, 6f8, 6f8, 7f8, 7f8, 5f8, 2i4, f8, 10f8, 10f8, 6f8, 500f8, 400f8, 200f8, 200f8, 360f8, 25f8, 10f8, 27f8, 9f8, 15f8, i4")
	sim_output_data = np.fromfile(f, dtype=dt, count=-1)
	f.close()


	lnL = np.zeros((len(sim_output_data),7))
	Param_Ep = np.zeros((len(sim_output_data), 6))

	for i in range(len(sim_output_data)):
		lnL[i] = extract_quantity(sim_output_data, 7, i)
		Param_Ep[i] = extract_quantity(sim_output_data, 5, i)


	ax.plot(Param_Ep[:,1], lnL[:,0]-np.max(lnL[:,0]), label='Total', marker='o', markersize=3)
	ax.plot(Param_Ep[:,1], lnL[:,3]-np.max(lnL[:,3]), label='Stern', marker='o', markersize=3)
	ax.plot(Param_Ep[:,1], lnL[:,5]-np.max(lnL[:,5]), label='EpGBM', marker='o', markersize=3)
	ax.plot(Param_Ep[:,1], 10*(lnL[:,6]-np.max(lnL[:,6])), label='eBAT6', marker='o', markersize=3)
	ax.set_ylabel(r'$\mathrm{ln(\mathcal{L}) - ln(\mathcal{L}_{max})}$', **font)
	ax.set_xlabel(r'$\mathrm{\sigma_{Ep}}$', **font)
	ax.legend()

	return

def plot_lnL_vs_alpha_amati(ax, run):
	binary_file = root_dir + run+'/reprise.dat'
	f = open(binary_file,'rb')
	dt = np.dtype("i4, 24i4, 10f8, 10f8, 6f8, 6f8, 7f8, 7f8, 5f8, 2i4, f8, 10f8, 10f8, 6f8, 500f8, 400f8, 200f8, 200f8, 360f8, 25f8, 10f8, 27f8, 9f8, 15f8, i4")
	sim_output_data = np.fromfile(f, dtype=dt, count=-1)
	f.close()


	lnL = np.zeros((len(sim_output_data),7))
	Param_Ep = np.zeros((len(sim_output_data), 6))

	for i in range(len(sim_output_data)):
		lnL[i] = extract_quantity(sim_output_data, 7, i)
		Param_Ep[i] = extract_quantity(sim_output_data, 5, i)


	ax.plot(Param_Ep[:,3], lnL[:,0]-np.max(lnL[:,0]), label='Total', marker='o', markersize=3)
	ax.plot(Param_Ep[:,3], lnL[:,3]-np.max(lnL[:,3]), label='Stern', marker='o', markersize=3)
	ax.plot(Param_Ep[:,3], lnL[:,5]-np.max(lnL[:,5]), label='EpGBM', marker='o', markersize=3)
	ax.plot(Param_Ep[:,3], 10*(lnL[:,6]-np.max(lnL[:,6])), label='eBAT6', marker='o', markersize=3)
	ax.set_ylabel(r'$\mathrm{ln(\mathcal{L}) - ln(\mathcal{L}_{max})}$', **font)
	ax.set_xlabel(r'$\mathrm{\alpha_{Amati}}$', **font)
	ax.legend()

	return


def plot_lnL_vs_zeta(ax, run):
	binary_file = root_dir + run+'/reprise.dat'
	f = open(binary_file,'rb')
	dt = np.dtype("i4, 24i4, 10f8, 10f8, 6f8, 6f8, 7f8, 7f8, 5f8, 2i4, f8, 10f8, 10f8, 6f8, 500f8, 400f8, 200f8, 200f8, 360f8, 25f8, 10f8, 27f8, 9f8, 15f8, i4")
	sim_output_data = np.fromfile(f, dtype=dt, count=-1)
	f.close()


	lnL = np.zeros((len(sim_output_data),7))
	Param_z = np.zeros((len(sim_output_data), 10))

	for i in range(len(sim_output_data)):
		lnL[i] = extract_quantity(sim_output_data, 7, i)
		Param_z[i] = extract_quantity(sim_output_data, 3, i)


	ax.plot(Param_z[:,7], lnL[:,0]-np.max(lnL[:,0]), label='Total', marker='o', markersize=3)
	ax.plot(Param_z[:,7], lnL[:,3]-np.max(lnL[:,3]), label='Stern', marker='o', markersize=3)
	ax.plot(Param_z[:,7], lnL[:,5]-np.max(lnL[:,5]), label='EpGBM', marker='o', markersize=3)
	ax.plot(Param_z[:,7], 10*(lnL[:,6]-np.max(lnL[:,6])), label='eBAT6', marker='o', markersize=3)
	lnLtot = lnL[:,5] + lnL[:,3] + lnL[:,6]
	ax.plot(Param_z[:,7], lnLtot-np.max(lnLtot), label='Total no weight', marker='o', markersize=3)
	ax.set_ylabel(r'$\mathrm{ln(\mathcal{L}) - ln(\mathcal{L}_{max})}$', **font)
	ax.set_xlabel(r'$\mathrm{\zeta}$', **font)
	ax.legend()

	return

def plot_param(param):
	fig = plt.figure(figsize=(12,10))
	ax = fig.add_subplot(111)

	run = '171215_lnl_vs_{}'.format(param)
	exec("plot_lnL_vs_{}(ax, run)".format(param))

	fig.savefig(outputdir+'lnL_vs_{}_newlnL.png'.format(param))
	ax.set_ylim(-6, 1)
	fig.savefig(outputdir+'lnL_vs_{}_Amati_newlnL_zoomed.png'.format(param))
	
	return


#plot_param(param)
fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(111)

run = '180104_PL_LogNorm_zeta_check_best'
plot_lnL_vs_zeta(ax, run)

ax.axvline(0.34)
ax.axvline(0.38)
ax.axvline(0.40)
#fig.savefig(outputdir+'lnL_vs_{}_newlnL.png'.format(param))
#ax.set_ylim(-6, 1)
#fig.savefig(outputdir+'lnL_vs_{}_Amati_newlnL_zoomed.png'.format(param))
	
plt.show()
