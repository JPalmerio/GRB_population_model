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
import subprocess

"""
This program launcher is used to initiate the GRB_population_MC_(14/16).exe
"""


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
        print("[INFO] Path '%s' created." % path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            print("[INFO] Path '%s' already exists." % path)


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    if v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def main():
    parser = argparse.ArgumentParser()
    r = subprocess.run(['sysctl', '-n', 'hw.ncpu'], stdout=subprocess.PIPE)
    ncmax = int(r.stdout)
    parser.add_argument("nb_procs", type=int, help=f"Number of processors used in code ({ncmax} max)")
    parser.add_argument("--nice", type=str2bool, nargs='?', const=True, default=True, help="Activate nice mode.")
    parser.add_argument("--out", type=str2bool, nargs='?', const=True, default=True, help="Save output to file.")
    args = parser.parse_args()

    og_input_directory = '../Input_para/'
    init_filename      = 'GRB_pop.init'
    init_filename_lum  = 'lum_param_search.init'
    init_filename_z    = 'redshift_param_search.init'
    init_filename_Ep   = 'Ep_param_search.init'
    init_filename_spec = 'spectrum_param_search.init'

    init_file_data = pf.read_column(og_input_directory+init_filename, 0, dtype=str)
    # Make sure output path exists
    path = '../Model_outputs/' + init_file_data[0]
    inputdir = 'Inputs/'
    make_sure_path_exists(path + inputdir)

    if int(init_file_data[3]) == 2:  # If post_processing mode enabled
        outputlog_file = 'outputlog_pp.txt'
        # Remove old reprise_pp file
        print("[INFO] You chose post-processing mode -> cleaning up old files.")
        os.system("rm  %s*post_proc.dat" % path)
        if os.path.isfile(path+'reprise_pp.dat'):
            os.remove(path+'reprise_pp.dat')
        # copy initfiles from '../Model_outputs/run/Inputs/' directory to '../Input_para/'
        os.remove(og_input_directory + init_filename     )
        os.remove(og_input_directory + init_filename_lum )
        os.remove(og_input_directory + init_filename_z   )
        os.remove(og_input_directory + init_filename_Ep  )
        os.remove(og_input_directory + init_filename_spec)

        print("[INFO] You chose post-processing mode : removed all input files from %s and copied them from %s" % (og_input_directory, path + inputdir))
        copyfile(path + inputdir + init_filename     , og_input_directory + init_filename     )
        copyfile(path + inputdir + init_filename_lum , og_input_directory + init_filename_lum )
        copyfile(path + inputdir + init_filename_z   , og_input_directory + init_filename_z   )
        copyfile(path + inputdir + init_filename_Ep  , og_input_directory + init_filename_Ep  )
        copyfile(path + inputdir + init_filename_spec, og_input_directory + init_filename_spec)

        # Change the value of run mode because copying from '../Model_outputs/run/Inputs/' changed that value to 1 (and wewant 2 for pp mode)
        f = open(og_input_directory + init_filename, 'r')
        lines = f.readlines()
        lines[23-1] = "     2\n"    # n is the line number you want to edit; subtract 1 as indexing of list starts from 0
        f.close()   # close the file and reopen in write mode to enable writing to file; you can also open in append mode and use "seek", but you will have some unwanted old data if the new data is shorter in length.

        f = open(og_input_directory + init_filename, 'w')
        f.writelines(lines)
        f.close()

    elif int(init_file_data[1]) == 1:  # If reprise mode
        outputlog_file = 'outputlog_reprise.txt'
        print("[INFO] You chose reprise mode.")
        inputdir = 'Inputs_reprise/'
        make_sure_path_exists(path + inputdir)
        # copy inputfiles from '../Input_para/' to '../Model_outputs/run/Inputs_reprise/' directory to keep a trace and for post processing
        copyfile(og_input_directory + init_filename     , path + inputdir + init_filename     )
        copyfile(og_input_directory + init_filename_lum , path + inputdir + init_filename_lum )
        copyfile(og_input_directory + init_filename_z   , path + inputdir + init_filename_z   )
        copyfile(og_input_directory + init_filename_Ep  , path + inputdir + init_filename_Ep  )
        copyfile(og_input_directory + init_filename_spec, path + inputdir + init_filename_spec)
    else:
        outputlog_file = 'outputlog.txt'

        # copy inputfiles from '../Input_para/' to '../Model_outputs/run/Inputs/' directory to keep a trace and for post processing
        copyfile(og_input_directory + init_filename     , path + inputdir + init_filename     )
        copyfile(og_input_directory + init_filename_lum , path + inputdir + init_filename_lum )
        copyfile(og_input_directory + init_filename_z   , path + inputdir + init_filename_z   )
        copyfile(og_input_directory + init_filename_Ep  , path + inputdir + init_filename_Ep  )
        copyfile(og_input_directory + init_filename_spec, path + inputdir + init_filename_spec)

    # prints the top 20 available computers with at least 4 procs and no users logged ranked by number of procs
    # os.system("taupe50 -t 20 -s proc_nb -q -n 4 -r")
    # possible computers to use : sideritis, pinot, picardan, romorantin, abouriou, cannes, colombard

    if args.out:
        print("[INFO] Saving output in %s%s" % (path, outputlog_file))
        if args.nice:
            print("[INFO] Executing nicely GRB_population_MC.exe on %d processors." % (args.nb_procs))
            os.system("nohup nice -n 19 mpiexec -n %d ./GRB_population_MC.exe > %s%s &" % (args.nb_procs, path, outputlog_file))
        else:
            print("[INFO] Executing GRB_population_MC.exe on %d processors." % (args.nb_procs))
            os.system("nohup mpiexec -n %d ./GRB_population_MC.exe > %s%s &" % (args.nb_procs, path, outputlog_file))
    else:
        print("[INFO] Sending output to terminal.")
        if args.nice:
            print("[INFO] Executing nicely GRB_population_MC.exe on %d processors." % (args.nb_procs))
            os.system("nice -n 19 mpiexec -n %d ./GRB_population_MC.exe " % (args.nb_procs))
        else:
            print("[INFO] Executing GRB_population_MC.exe on %d processors." % (args.nb_procs))
            os.system("mpiexec -n %d ./GRB_population_MC.exe " % (args.nb_procs))
    return


if __name__ == "__main__":
    main()
