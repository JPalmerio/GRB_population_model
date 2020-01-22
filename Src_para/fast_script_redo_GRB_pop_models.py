# To allow the importation of plotting function module anywhere
import sys
import platform
if platform.system() == 'Linux':
	sys.path.insert(0,'/nethome/palmerio/Dropbox/Plotting_GUI/Src')
elif platform.system() == 'Darwin': 
	sys.path.insert(0,'/Users/palmerio/Dropbox/Plotting_GUI/Src')
import os
import errno
import plotting_functions as pf
from shutil import copyfile
import argparse
import fileinput

mod = ['A', 'LN']
k = ['0', '2']
sig = '6.5'
core_name = '180511_SVOM_'

path = '/home/versailles1NS/palmerio/Model_outputs/'
path_for_sim = '/nethome/palmerio/1ere_annee/Frederic/GRB_population_code/Input_para/'
path_for_output = '/nethome/palmerio/1ere_annee/Frederic/GRB_population_code/Model_outputs/'

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
        print "[INFO] Path '%s' created." %path
    except OSError as exception:
        if exception.errno != errno.EEXIST:
           raise
        else :
           print "[INFO] Path '%s' already exists." %path
    return

def create_paths():
	path_list = []
	for i in range(len(mod)):
		for j in range(len(k)):
			path_list.append(path+core_name+mod[i]+'_k'+k[j]+'_sig'+sig)
			print path_list[-1].split('/')[-1]
	return path_list

path_list = create_paths()
#for i in range(len(path_list)):
i = 2
make_sure_path_exists(path_for_output + path_list[i].split('/')[-1])
outputlog_file = path_for_output + path_list[i].split('/')[-1]+'/outputlog.txt'
print 'Recomputing model : {} into {}'.format(path_list[i], path_for_sim + 'GRB_pop.init')
copyfile(path_list[i]+'/Inputs/GRB_pop.init' ,  path_for_sim + 'GRB_pop.init'   )
os.system("python code_launcher.py 16 1")





