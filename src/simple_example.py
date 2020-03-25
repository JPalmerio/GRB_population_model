from ECLAIRs import init_ECLAIRs
from cosmology import init_cosmology
from GRB_population import GRBPopulation
import io_grb_pop as io
import miscellaneous as msc
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
incl_samples, incl_instruments, incl_constraints = msc.create_config(config,
                                                                     samples,
                                                                     instruments,
                                                                     obs_constraints)
# Generate the GRB population
np.random.seed(0)
GRB_pop = GRBPopulation(Nb_GRBs=int(float(config['Nb_GRBs'])),
                        output_dir=paths_to_dir['output'])
GRB_pop.draw_GRB_properties(cosmo=cosmo, params=params, run_mode='debug', savefig=True)
GRB_pop.calc_peak_photon_flux(incl_instruments=incl_instruments)
GRB_pop.calc_det_prob(incl_samples=incl_samples)
GRB_pop.create_mock_constraints(obs_constraints=incl_constraints)
GRB_pop.compare_to_observational_constraints(obs_constraints=incl_constraints, method='chi2')

log.info(GRB_pop.summary())
