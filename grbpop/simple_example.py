from ECLAIRs import init_ECLAIRs
from cosmology import init_cosmology, create_cosmology
from GRB_population import GRBPopulation, create_GRB_population_from
import io_grb_pop as io
import numpy as np
import logging
import sys

log = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                    format='%(asctime)s.%(msecs)03d [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
logging.getLogger('matplotlib').setLevel(logging.WARNING)


# Define the paths used by the code
paths_to_dir, paths_to_files = io.generate_paths(conf_fname='config_simple_example.yml',
                                                 param_fname='parameters_simple_example.yml',
                                                 init_dir=None)

# Read the input files
config, params, instruments, samples, obs_constraints = io.read_init_files(paths_to_files)
io.create_output_dir(paths_to_dir=paths_to_dir, dir_name=config['output_dir'], overwrite=True)

# If you want to make predictions for ECLAIRs you need the ECLAIRs properties
ECLAIRs_prop = None
# ECLAIRs_prop = init_ECLAIRs(ECLAIRs_dir=paths_to_dir['ECLAIRs'],
#                             ECLAIRs_config=instruments['ECLAIRs'])
# samples['ECLAIRs']['pflx_min'] = ECLAIRs_prop['bkg_total']
# You must do this before calling io.create_config

# Code calculates which samples, instruments, and constraints to include
incl_samples, incl_instruments, incl_constraints = io.create_config(config,
                                                                    samples,
                                                                    instruments,
                                                                    obs_constraints)

# Initialize the cosmology
cosmo = init_cosmology(paths_to_dir['cosmo'])
# If you want to create your own cosmology use:
# cosmo = create_cosmology(OmegaM=0.3, OmegaL=0.7, h=0.7)

# Generate the GRB population
np.random.seed(0)
gp = create_GRB_population_from(Nb_GRBs=config['Nb_GRBs'],
                                cosmo=cosmo,
                                params=params,
                                incl_samples=incl_samples,
                                incl_instruments=incl_instruments,
                                incl_constraints=incl_constraints,
                                ECLAIRs_prop=ECLAIRs_prop,
                                output_dir=paths_to_dir['output'],
                                run_mode='debug',
                                savefig=True)
# This function is equivalent to the following:
# gp = GRBPopulation(Nb_GRBs=config['Nb_GRBs'], output_dir=paths_to_dir['output'])
# gp.draw_GRB_properties(cosmo=cosmo, params=params, run_mode='debug', savefig=True)
# gp.calculate_quantities(instruments=incl_instruments, samples=incl_samples)
# gp.create_mock_constraints(constraints=incl_constraints)
# gp.compare_to_observational_constraints(constraints=incl_constraints)
# gp.normalize_to_Stern()
# print(gp.summary())
