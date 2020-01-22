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

plt.style.use('ggplot')
homedir = '/nethome/palmerio/1ere_annee/Frederic/'
text_size = 22
font = {'fontname':'Serif', 'size':text_size}

#filename = homedir + 'catalogs/Swift_cat/Swift_redshift_cat.txt'
#outputfile = 'Swift_temp_save_table.txt'
def clean_Swift_redshift_database(filename, outputfile, name_row=0, redshift_row=4):
	"""
		A function that scans through a filename database and asks the user to input the correct redshift value and rewrites a clean version in a new file.
	"""
	data = pd.read_csv(filename, sep='\t', skiprows=1)
	
	name = data.values.transpose()[name_row]
	table_redshift = data.values.transpose()[redshift_row] 
	
	name_final = []
	redshift_final = []
	
	init = int(input('Please enter the row at which to start (0-indexed) :\n'))
	
	save_file = homedir+'catalogs/Swift_cat/' + outputfile 
	print "The GRB list is %d lines long." %len(table_redshift)
	print "Please enter the final redshift for each burst : "
	for i in range(init,len(table_redshift)):
		if type(table_redshift[i]) == str:
			print i," at ",(int(100*i/len(table_redshift)))," %", name[i],' : ', table_redshift[i]
		#try:
			z = float(input())
			redshift_final.append(z)
			name_final.append(name[i])
		#except NameError, SyntaxError, ValueError:
		#	pass
			f = open(save_file, 'a')
			f.write("%s\t %6.4lf\n" %(name_final[-1], redshift_final[-1]))
			f.close()
	return


# Swift data
filename = homedir + 'catalogs/Swift_cat/Swift_redshift_cleaned_table.txt'
data = pd.read_csv(filename, delimiter='\t',skiprows=1)
Swift_name = data.values.transpose()[0]
Swift_redshift = data.values.transpose()[1]

# eBAT6 data
filename = homedir + 'catalogs/BAT6_cat/eBAT6_cat.txt'
eBAT6_name     = pf.read_column(filename, 0, dtype=str, splitter='\t|', stripper='|') 
eBAT6_redshift = pf.read_column(filename, 1, splitter='\t|', stripper='|') 
eBAT6_redshift_masked = eBAT6_redshift[np.isfinite(eBAT6_redshift)]


fig = plt.figure(tight_layout=True, figsize=(10,8))
ax  = fig.add_subplot(111)

pf.ks_and_plot(ax, Swift_redshift, eBAT6_redshift_masked, labels=['Swift','eBAT6'],lw=2)

#Swift_redshift_for_ECDF, ECDF_Swift = pf.unbinned_empirical_cdf(Swift_redshift)
#eBAT6_redshift_for_ECDF, ECDF_eBAT6 = pf.unbinned_empirical_cdf(eBAT6_redshift_masked)
# ax.plot(Swift_redshift_for_ECDF, ECDF_Swift,drawstyle='steps-post', lw=2, label='Swift')
# ax.plot(eBAT6_redshift_for_ECDF, ECDF_eBAT6,drawstyle='steps-post', lw=2, label='eBAT6')

ax.legend(loc='best')
ax.set_xlabel(r'Redshift (z)', **font)
ax.set_ylabel(r'ECDF', **font)
ax.set_title('Empirical redshift cumulative distribution function\n of eBAT6 sample and Swift catalog', **font)

plt.show()