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
#from astroML.plotting import hist
import seaborn as sns
import scipy as sp
import datetime
import pandas as pd
from matplotlib.dates import date2num
from matplotlib.dates import num2date

#from pyqt_fit import kde, kde_methods
plt.style.use('default')

root_dir = '/nethome/palmerio/1ere_annee/Frederic/'
GRB_pop_dir = 'GRB_population_code/Model_outputs/'
text_size = 22
font = {'fontname':'Serif', 'size':text_size}

def replace_nans(obs_redshift, z_unknown):
	obs_redshift_masked = np.zeros(len(obs_redshift))
	for i in range(len(obs_redshift)):
		if np.isfinite(obs_redshift[i]):
			obs_redshift_masked[i] = obs_redshift[i]
		else:
			obs_redshift_masked[i] = z_unknown
	return obs_redshift_masked

def ks_and_plot(ax, sample_1, sample_2, xy=(0.5, 0.3), colors=[None,None], labels=[None,None], title=None, **kwargs):

	ks_D_stat, p_value = sp.stats.ks_2samp(sample_1, sample_2)
	unbinned_sample_1, ECDF_1 = pf.unbinned_empirical_cdf(sample_1)
	unbinned_sample_2, ECDF_2 = pf.unbinned_empirical_cdf(sample_2)
	if colors[0] is None :
		#ax.hist(sample_1, bins=500,cumulative=True, histtype='step', normed=True,label=labels[0], **kwargs)
		ax.plot(unbinned_sample_1, ECDF_1 ,drawstyle='steps-post',label=labels[0], **kwargs)
	else :
		#ax.hist(sample_1, bins=500,cumulative=True, histtype='step', normed=True,label=labels[0],color=colors[0], **kwargs)
		ax.plot(unbinned_sample_1, ECDF_1 ,drawstyle='steps-post',label=labels[0],color=colors[0], **kwargs)

	if colors[1] is None :
		#ax.hist(sample_2, bins=500,cumulative=True, histtype='step', normed=True,label=labels[1], **kwargs)
		ax.plot(unbinned_sample_2, ECDF_2 ,drawstyle='steps-post',label=labels[1], **kwargs)

	else :
		#ax.hist(sample_2, bins=500,cumulative=True, histtype='step', normed=True,label=labels[1],color=colors[1], **kwargs)
		ax.plot(unbinned_sample_2, ECDF_2 ,drawstyle='steps-post',label=labels[1],color=colors[1], **kwargs)

	ax.annotate("KS D = %.4lf\np_value = %.4lf"%(ks_D_stat, p_value), xy=xy, xycoords='axes fraction')
	ax.legend(loc='best')
	ax.set_xlabel(r'Redshift z')
	ax.set_ylabel(r'Cumulative Distribution Function')
	ax.set_xlim(0,max(max(sample_1), max(sample_2)))
	if title is None:
		ax.set_title(r'Redshift distribution for the extended BAT6 sample')
	else :
		ax.set_title(title)

	return

def easy_plot(ax, x, y, colormap=None, x_is_log=False, y_is_log=False, capsize=0, errorcolor='k', **kwargs):
	"""
	Helper function to easily plot two variables with incomplete data sets (uses masks) and errors, or limits.
	Assumes the data is in the form of the output of read_column.
	i.e. : x[0] = data (float)
		   x[1] = error plus (float)
		   x[2] = error minus (float)
		   x[3] = upper limit (bool)
		   x[4] = lower limit (bool)
	Colormap input is expected as such :
		Colormap[0] is an array of values assumed to be of same length as x and y.
		Colormap[1] is one value (lower value for colorbar scale)
		Colormap[2] is one value (upper value for colorbar scale)
		Colormap[3] is the argument to put in cmap (ex: 'YlOrRd_r')
	ax : axes object on which to plot.
	Returns scatter plot artist (used to create colorbar afterward). 
	"""

	# Create masks
	x_mask = np.isfinite(x[0])
	y_mask = np.isfinite(y[0])

	# filters data
	if x_is_log :
		x = log_to_lin(x)
		ax.set_xscale('log')

	x_to_plot = x[0][x_mask & y_mask]
	xerr = np.asarray([x[2][x_mask & y_mask], x[1][x_mask & y_mask]])
	xuplims = x[3][x_mask & y_mask]
	xlolims = x[4][x_mask & y_mask]
	for i in range(len(x_to_plot)):
		if xuplims[i]:
			xerr[0][i] = 0.5 * x_to_plot[i] # lower error for upper limit (arrow pointing down)
		if xlolims[i] :
			xerr[1][i] = 1.5 * x_to_plot[i] # upper error for lower limit (arrow pointing up)
		# The coefficients were chosen to give sizeable arrows in log scale

	if y_is_log:
		y = log_to_lin(y)
		ax.set_yscale('log')

	y_to_plot = y[0][x_mask & y_mask]
	yerr = np.asarray([y[2][x_mask & y_mask], y[1][x_mask & y_mask]])
	yuplims = y[3][x_mask & y_mask]
	ylolims = y[4][x_mask & y_mask]
	for i in range(len(y_to_plot)):
		if yuplims[i]:
			yerr[0][i] = 0.5 * y_to_plot[i] # lower error for upper limit (arrow pointing down)
		if ylolims[i] :
			yerr[1][i] = 1.5 * y_to_plot[i] # upper error for lower limit (arrow pointing up)
		# The coefficients were chosen to give sizeable arrows in log scale


	# Plotting
	ax.errorbar(x_to_plot, y_to_plot, xerr=xerr, yerr=yerr, xuplims=xuplims, xlolims=xlolims, uplims=yuplims, lolims=ylolims, capsize=capsize, color=errorcolor, marker=None, fmt='.', zorder=0)
	if colormap != None :
		norm = matplotlib.colors.Normalize(vmin=colormap[1], vmax=colormap[2]) # limits to the colorbar if colormap is used
		scatterplot = ax.scatter(x_to_plot, y_to_plot, c=colormap[0][x_mask & y_mask], cmap=colormap[3], norm=norm, **kwargs)
	else:
		scatterplot = ax.scatter(x_to_plot, y_to_plot, **kwargs)

	return scatterplot


# Create original BAT6
filename_og = root_dir + 'catalogs/BAT6_cat/BAT6_2012.txt'
name_og  = pf.read_column(filename_og,  0, dtype=str, splitter = '\t|')
redshift_og  = pf.read_column(filename_og,  1, dtype=float, splitter = '\t|')
redshift_og_mask = np.isfinite(redshift_og)
obs_redshift_og_masked = np.zeros(len(redshift_og))


# Extended BAT6 : eBAT6
file_eBAT6_obs = root_dir + 'catalogs/BAT6_cat/eBAT6_cat.txt'
obs_name     = pf.read_column(file_eBAT6_obs, 0, stripper='|', splitter = '|', dtype=str)
obs_redshift = pf.read_column(file_eBAT6_obs, 1, stripper='|', splitter = '|')
obs_redshift2 = pf.read_data(file_eBAT6_obs, 1, stripper='|', splitter = '|')
obs_redshift_mask = np.isfinite(obs_redshift)
obs_redshift_masked = np.zeros(len(obs_redshift))
obs_redshift_masked = obs_redshift[obs_redshift_mask]

file_Swift_obs = root_dir + 'catalogs/Swift_cat/Swift_pflx_cat.txt'
obs_name_S2     = pf.read_column(file_Swift_obs, 0, dtype=str)
obs_name_S = []
for i in range(len(obs_name_S2)):
	if int(obs_name_S2[i][:2]) <= 14 :
		obs_name_S.append(obs_name_S2[i])
obs_name_S = np.asarray(obs_name_S).astype(str)

file_Swift_obsb = root_dir + 'catalogs/Swift_cat/Swift_cat.txt'
obs_name_S2b     = pf.read_column(file_Swift_obsb, 0, dtype=str)
obs_name_Sb = []
for i in range(len(obs_name_S2b)):
	if int(obs_name_S2b[i][:2]) <= 14 :
		obs_name_Sb.append(obs_name_S2b[i])
obs_name_Sb = np.asarray(obs_name_Sb).astype(str)



# header = pf.read_overhead(file_eBAT6_obs, stripper='|', splitter = '|')
# for i in range(len(header)):
# 	print i, header[i]

# bins = np.linspace(0,6,16)
# bins_mid = 0.5*(bins[1:]+bins[:-1])
# hist_eBAT6, bins_BAT6 = np.histogram(obs_redshift_masked, bins=bins)
# err = pf.poisson_errors(hist_eBAT6)

def convert_GRBname_to_days(date_string):
	year = int(date_string[:2]) + 2000
	month = int(date_string[2:4])
	day = int(date_string[4:6])
	float_days = date2num(datetime.datetime(year, month, day))

	return float_days


def GRB_date_histogram(GRB_names, bins=None, ax=None, show_plot=False, **kwargs):

	GRB_names_in_days = np.zeros(len(GRB_names))
	for i in range(len(GRB_names)):
		GRB_names_in_days[i] = convert_GRBname_to_days(GRB_names[i])

	if bins is None:
		bins = np.linspace(convert_GRBname_to_days('050101'),convert_GRBname_to_days('150101'), 11)

	hist, bins = np.histogram(GRB_names_in_days, bins=bins)
	if show_plot:
		if ax is None:
			ax = plt.gca()
		ax.hist(GRB_names_in_days, bins=bins, **kwargs)

	return hist, bins



fig = plt.figure(tight_layout=True, figsize=(10,8))
ax = fig.add_subplot(111)
ax.set_title('Fraction of Swift LGRBs in eBAT6 sample', **font)
hist_S, bins = GRB_date_histogram(obs_name_S, ax=ax, show_plot=True, label='All Swift', linewidth=1, edgecolor='k')
hist_Sb, bins = GRB_date_histogram(obs_name_Sb, ax=ax, show_plot=True, label='Swift w/ redshift', linewidth=1, edgecolor='k')
hist, bins = GRB_date_histogram(obs_name, ax=ax, show_plot=True, label='eBAT6', linewidth=1, edgecolor='k')

eBAT6_percentage = 100*hist.astype(float)/hist_S.astype(float)
dates_x = 0.5*(bins[:-1]+bins[1:])
for i in range(len(dates_x)):
	ax.text(dates_x[i], eBAT6_percentage[i], '{:.1f}%'.format(eBAT6_percentage[i]), ha='center', size=15)

ax.set_xticks(dates_x)
ax.set_xticklabels([2005+i for i in range(10)])
# bins_in_date = num2date(bins)
# bins_in_date_label = []
# for i in range(len(bins_in_date)):
# 	bins_in_date_label.append(bins_in_date[i].year)
# ax.set_xticklabels(bins_in_date_label)

for tick in ax.get_xticklabels():
    tick.set_rotation(45)
ax.legend()

ax.set_xlabel('Year', **font)
ax.set_ylabel('Number of GRBs', **font)
plt.show()

raise SystemExit

# names_dt = []
# for i in range(len(obs_name_S)):
# 	names_dt.append(convert_to_datetime_object(obs_name_S[i]))
# df = pd.DataFrame({'date':names_dt})
# df["date"] = df["date"].astype("datetime64")
# df.groupby(df['date'].dt.year).count().plot(ax=plt.gca(), color='g', alpha=0.7,kind='bar', label='All Swift')

# names_dt = []
# for i in range(len(obs_name_Sb)):
# 	names_dt.append(convert_to_datetime_object(obs_name_Sb[i]))
# df = pd.DataFrame({'date':names_dt})
# df["date"] = df["date"].astype("datetime64")
# df.groupby(df['date'].dt.year).count().plot(ax=plt.gca(), color='darkorange', alpha=0.9,kind='bar', label='Swift w/ redshift')


# names_dt = []
# for i in range(len(obs_name)):
# 	names_dt.append(convert_to_datetime_object(obs_name[i]))
# df = pd.DataFrame({'date':names_dt})
# df["date"] = df["date"].astype("datetime64")
# df.groupby(df['date'].dt.year).count().plot(ax=plt.gca(),kind='bar', label='eBAT6')


plt.show()

raise SystemExit

write_catalog = False
if write_catalog :
	f = open(root_dir +'GRB_population_code/observational_constraints/eBAT6_constraint.txt', 'w')
	f.write("# eBAT6 redshift distribution constraint from Pescalli et al. 2016\n")
	f.write("# redshift bin | counts | poisson error\n")
	f.write("# Number of GRBs in sample : %d\n" %sum(hist_eBAT6))
	for i in range(len(hist_eBAT6)):
		f.write("%e \t %d \t %d\n" %(bins_BAT6[i], hist_eBAT6[i], err[0][i]))
		f.write("%e \t %d \t %d\n" %(bins_BAT6[i+1], hist_eBAT6[i], err[0][i]))
	f.close()

plot_redshift = False
if plot_redshift : 
	# Plot stuff
	fig = plt.figure(tight_layout=True, figsize=(16,12))
	ax1 = fig.add_subplot(221)
	ax2 = fig.add_subplot(222)
	ax3 = fig.add_subplot(223)
	ax4 = fig.add_subplot(224)

	#hist(obs_redshift_masked, bins='blocks',   label ='Bayesian blocks', alpha=0.95,lw=2,  histtype='step', normed=True)
	##True = kde.KDE1D(obs_redshift_masked)
	#hist(obs_redshift_masked, bins='knuth',    label ='Knuth',           alpha=0.95,lw=2,  histtype='step', normed=True)
	#hist(obs_redshift_masked, bins='scott',    label ='Scott',           alpha=0.95,lw=2,  histtype='step', normed=True)
	#hist(obs_redshift_masked, bins='freedman', label ='Freedman',        alpha=0.95,lw=2,  histtype='step', normed=True)
	obs_redshift_masked = obs_redshift[obs_redshift_mask]
	obs_redshift_og_masked = redshift_og[redshift_og_mask]
	
	obs_redshift_og_masked1 = replace_nans(redshift_og, 0.)
	obs_redshift_og_masked2 = replace_nans(redshift_og, 8.)
	ks_and_plot(ax1, obs_redshift_og_masked1, obs_redshift_og_masked, colors=['darkorange','k'], xy=(0.3,0.3), labels=['z=0 OG', None], alpha=0.95, lw=2)
	ks_and_plot(ax1, obs_redshift_og_masked2, obs_redshift_og_masked, colors=['purple','k'], xy=(0.3,0.5), labels=['z=8 OG','OG Bat6'], alpha=0.95, lw=2, title='Original BAT6 z distr. assuming missing data is at z=0 or z=8')
	
	ks_and_plot(ax2, obs_redshift_masked, obs_redshift_og_masked ,colors=['darkorange','purple'], labels=['eBat6','OG Bat6'], alpha=0.95, lw=2, title='Original and extended BAT6 z distr. compared')
	
	z_unknown = 0.
	obs_redshift_masked = replace_nans(obs_redshift, z_unknown)
	obs_redshift_og_masked = replace_nans(redshift_og, z_unknown)
	ks_and_plot(ax3, obs_redshift_masked, obs_redshift_og_masked ,colors=['darkorange','purple'], labels=['z=0 ext','z=0 og'], alpha=0.95, lw=2, title='Original and extended BAT6 z distr. assuming missing data is at z=0')
	
	z_unknown = 8.
	obs_redshift_masked = replace_nans(obs_redshift, z_unknown)
	obs_redshift_og_masked = replace_nans(redshift_og, z_unknown)
	ks_and_plot(ax4, obs_redshift_masked, obs_redshift_og_masked ,colors=['darkorange','purple'], labels=['z=8 ext','z=8 og'], alpha=0.95, lw=2, title='Original and extended BAT6 z distr. assuming missing data is at z=8')
	
	
	#hist(obs_redshift_masked, bins= 10,        label ='10 bins',         alpha=0.95,lw=2,  histtype='step', normed=True)
	#sns.kdeplot(obs_redshift_masked, ax=ax, label='KDE', lw=2, color='k')
	#sns.rugplot(obs_redshift_masked, ax=ax)



plt.show()








