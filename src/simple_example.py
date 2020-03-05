from ECLAIRs import init_ECLAIRs
from cosmology import init_cosmology
from GRB_population import GRBPopulation
import io_grb_pop as io
import miscellaneous as msc
import observational_constraints as obs
import numpy as np
import logging
import sys

log = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                    format='%(asctime)s.%(msecs)03d [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
logging.getLogger('matplotlib').setLevel(logging.WARNING)


# Define the paths used by the code
paths_to_dir, paths_to_files = io.generate_paths(init_dir=None)
# Just change the paths to the config and param files
paths_to_files['config'] = paths_to_dir['init']/'config_simple_example.yml'
paths_to_files['param'] = paths_to_dir['init']/'parameters_simple_example.yml'

# Read the input files
config, params, instruments, samples, obs_constraints = io.read_init_files(paths_to_files)
io.create_output_dir(paths_to_dir=paths_to_dir, dir_name=config['output_dir'], overwrite=True)

# Initialize the cosmology, ECLAIRs properties and the observational
# constraints by reading the appropriate files
cosmo = init_cosmology(paths_to_dir['cosmo'])
ECLAIRs_prop = init_ECLAIRs(ECLAIRs_dir=paths_to_dir['ECLAIRs'],
                            ECLAIRs_config=instruments['ECLAIRs'])
samples['ECLAIRs']['pflx_min'] = ECLAIRs_prop['bkg_total']

# Some final initializations
obs_constraints = obs.load_observational_constraints(obs_constraints)
incl_samples = msc.included_samples(config['samples'], samples)
incl_instruments = msc.included_instruments(incl_samples, instruments)

# Generate the GRB population
np.random.seed(0)
GRB_population = GRBPopulation(Nb_GRBs=int(float(config['Nb_GRBs'])),
                               output_dir=paths_to_dir['output'])
GRB_population.draw_GRB_properties_for_MCMC(cosmo=cosmo,
                                            params=params,
                                            incl_instruments=incl_instruments)
GRB_population.create_mock_constraint(obs_constraints=obs_constraints)
GRB_population.compare_to_observational_constraints(obs_constraints=obs_constraints, method='chi2')

