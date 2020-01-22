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
from astroML.plotting import hist
#import seaborn as sns
import scipy as sp



root_dir = '/nethome/palmerio/1ere_annee/Frederic/'
filename = root_dir + 'catalogs/BAT6_cat/BAT6_2012.txt'

name = pf.read_column(filename, 0, dtype=str, splitter = '\t|')
redshift = pf.read_column(filename, 1, dtype=float, splitter = '\t|')
RA = pf.read_column(filename, 2, dtype=str, splitter = '\t|')
dec = pf.read_column(filename, 3, dtype=str, splitter = '\t|')

j=0
redshift_lim = 2
print "Original BAT6 sample (Salvaterra et al. 2012) northern (dec >= 30 degrees) GRBs :"
print "Redshift cut applied : %.3lf" %redshift_lim
print "no \t name   \t redshift \t RA     \t dec"
for i in range(len(name)):
	if float(dec[i][:3]) >= 30 :
		if redshift[i] <= redshift_lim:
			j += 1 
			print "%2d \t %s \t %.4lf \t %s \t %s" %(j, name[i], redshift[i], RA[i], dec[i])




file_eBAT6_obs = root_dir + 'catalogs/BAT6_cat/BAT6_cat.txt'
eBAT6_name = pf.read_column(file_eBAT6_obs, 0, dtype=str, stripper='|', splitter = '\t|')
eBAT6_redshift = pf.read_column(file_eBAT6_obs, 1, stripper='|', splitter = '|')
eBAT6_redshift_mask = np.isfinite(eBAT6_redshift)
eBAT6_redshift_masked = np.zeros(len(eBAT6_redshift))


file_swift_total = root_dir + 'catalogs/Swift_cat/Swift_BAT6_RA_dec_cat.txt'
#overhead = pf.read_overhead(file_swift_total, splitter='|', stripper='|')
#for i in range(len(overhead)):
#		print i, overhead[i]
for i in range(len(eBAT6_name)):
	eBAT6_name[i] = eBAT6_name[i].strip()

Swift_name = pf.read_column(file_swift_total, 0, dtype=str, splitter='|', stripper='|')
Swift_RA = pf.read_column(file_swift_total, 1, dtype=str, splitter='|', stripper='|')
Swift_dec = pf.read_column(file_swift_total, 2, dtype=str, splitter='|', stripper='|')
for i in range(len(Swift_name)):
	Swift_name[i] = Swift_name[i][4:].strip()
	 #Swift_dec[i] = Swift_dec[i][:3]

eBAT6_unobserved_sample = []

for j in range(len(eBAT6_name)):
	in_old_flag = False
	for k in range(len(name)):
		if eBAT6_name[j] == name[k]:
			in_old_flag = True
	if in_old_flag:
		pass	
	else :
		for i in range(len(Swift_name)):
			if Swift_name[i] == eBAT6_name[j]:
				if float(Swift_dec[i][:3]) <= 30:
					#print Swift_name[i], eBAT6_name[j], eBAT6_redshift[j], Swift_RA[i], Swift_dec[i]
					eBAT6_unobserved_sample.append([Swift_name[i], eBAT6_name[j], eBAT6_redshift[j], Swift_RA[i], Swift_dec[i]])

eBAT6_unobserved_sample = np.asarray(eBAT6_unobserved_sample).transpose()

print ""
print "Extended BAT6 sample (Pescalli et al. 2016) southern (dec <= 30 degrees) GRBs :"
#print "Redshift cut applied : %.3lf" %redshift_lim
print "no \t name   \t redshift \t RA     \t dec"
j=0
for i in range(len(eBAT6_unobserved_sample)):
	if float(eBAT6_unobserved_sample[4][i][:3]) <= 30 :
		#if eBAT6_unobserved_sample[2][i] <= redshift_lim:
		j += 1 
		print "%2d \t %s \t %s \t %s \t %s" %(j, eBAT6_unobserved_sample[0][i], eBAT6_unobserved_sample[2][i], eBAT6_unobserved_sample[3][i], eBAT6_unobserved_sample[4][i])
print "Note : there are 24 bursts missing (from 2013 on)"

