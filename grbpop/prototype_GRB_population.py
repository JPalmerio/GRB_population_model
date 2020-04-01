import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import logging
import time
import sys
import pandas as pd
import physics as ph
import io_grb_pop as io
import miscellaneous as msc
from GRB_population import GRBPopulation
from cosmology import init_cosmology
from ECLAIRs import init_ECLAIRs

matplotlib.rc('text', usetex=False)

# Set up logging
log = logging.getLogger(__file__.rsplit('/')[0])
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                    format='%(asctime)s.%(msecs)03d [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
logging.getLogger('matplotlib').setLevel(logging.WARNING)


""" This is a lightweight version of the fortran code to generate a GRB population """

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description="GRB population code")
    argparser.add_argument('run_mode', type=int, help='Run mode: 0 for one_run, 1 for MCMC')
    argparser.add_argument('--debug', type=msc.str2bool, nargs='?',
                           const=True, help='Boolean, activates debug mode')
    args = argparser.parse_args()

    print(args.debug)

    if args.debug:
        run_mode = 'debug'
        verbose = True
        log.setLevel(logging.DEBUG)
    else:
        run_mode = None
        verbose = False
        log.setLevel(logging.INFO)

    paths_to_dir, paths_to_files = io.generate_paths()
    config, params, instruments, samples = io.read_init_files(paths_to_files)
    cosmo = init_cosmology(paths_to_dir['cosmo'])
    ECLAIRs_prop = init_ECLAIRs(ECLAIRs_dir=paths_to_dir['ECLAIRs'],
                                ECLAIRs_config=instruments['ECLAIRs'])
    samples['ECLAIRs']['pflx_min'] = ECLAIRs_prop['bkg_total']

    Nb_GRBs = int(float(config['Nb_GRBs']))
    log.info(f"Number of GRBs simulated: {Nb_GRBs:.2e}")

    incl_samples = msc.included_samples(config['samples'], samples)
    incl_instruments = msc.included_instruments(incl_samples, instruments)

    if args.run_mode == 0:
        io.create_output_dir(paths_to_dir, dir_name=config['output_dir'],
                             run_mode=run_mode)
        np.random.seed(0)
        tstart = time.time()
        GRB_population = GRBPopulation(Nb_GRBs, output_dir=paths_to_dir['output'])
        GRB_prop = GRB_population.draw_GRB_properties(cosmo=cosmo, params=params)

        ph.calc_peak_photon_flux(GRB_prop, incl_instruments, ECLAIRs_prop)
        ph.calc_peak_energy_flux(GRB_prop, incl_instruments, ECLAIRs_prop)
        ph.calc_photon_fluence(GRB_prop, incl_instruments)
        ph.calc_energy_fluence(GRB_prop, incl_instruments)
        ph.calc_det_prob(GRB_prop, incl_samples, **ECLAIRs_prop)

        df = pd.DataFrame(GRB_prop)

        if config['save_all_GRBs']:
            GRB_population.save_all_GRBs(paths_to_dir['output'])
        tend = time.time()
        log.info(f"Execution time {tend-tstart:.3f}s")

    elif args.run_mode == 1:
        io.create_output_dir(paths_to_dir, dir_name=config['output_dir'],
                             run_mode=run_mode)
        np.random.seed(0)
        tstart = time.time()
        raise NotImplementedError("Sorry I'm lazy")
        tend = time.time()
        log.info(f"MC routine execution time {tend-tstart:.3f}s")
        exit()
    else:
        log.error('Invalid run mode')
    # end of code
    if args.debug:
        plt.show()
