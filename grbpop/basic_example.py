from cosmology import init_cosmology
from GRB_population import create_GRB_population_from
import io_grb_pop as io
import numpy as np
import logging
import sys

log = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(asctime)s.%(msecs)03d [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
logging.getLogger('matplotlib').setLevel(logging.WARNING)


# Define the paths used by the code
paths_to_dir, paths_to_files = io.generate_paths(conf_fname='config_basic_example.yml',
                                                 param_fname='parameters_simple_example.yml',
                                                 init_dir=None)

# Read the input files
config, params, instruments, samples, obs_constraints = io.read_init_files(paths_to_files)

# Code calculates which samples, instruments, and constraints to include
incl_samples, incl_instruments, incl_constraints = io.create_config(config,
                                                                    samples,
                                                                    instruments,
                                                                    obs_constraints)

# Initialize the cosmology
cosmo = init_cosmology(paths_to_dir['cosmo'])

# Generate the GRB population
np.random.seed(0)
gp = create_GRB_population_from(Nb_GRBs=config['Nb_GRBs'],
                                cosmo=cosmo,
                                params=params,
                                incl_samples=incl_samples,
                                incl_instruments=incl_instruments,
                                incl_constraints=incl_constraints,
                                output_dir=paths_to_dir['output'])
