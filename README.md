# GRB Population model

`grbpop` allows one to create a population of Gamma-Ray Bursts (GRBs).
Using Monte Carlo sampling from input distributions, one can then compare the predictions for various observational constraints.

## Basic usage
The arborescence is important for this software.
Make sure you have the following structure:
```
.
├── data
│   ├── ECLAIRs (optional)
│   └── cosmology (optional)
├── init
├── model_outputs
├── observational_constraints
└── grbpop
```
The `data/ECLAIRs` directory is optional; you only need it if you want to make prediction for ECLAIRs.
The `data/cosmology` directory is also optional; you can create your own cosmology if you wish (see this).

For this basic example you will want to start with:
```
import io_grb_pop as io
from GRB_population import GRBPopulation, create_GRB_population_from
import logging
log = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                    format='%(asctime)s.%(msecs)03d [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
```
You can change the logging level to `INFO` if you wish.

### Setting up the paths
The first thing to do is set up the paths to the various files so that the code knows where to look for what:
```
paths_to_dir, paths_to_files = io.generate_paths()
```
If you did not follow the directory structure presented above, you can provide a `init_dir` argument to tell the code where to look for the initialization files.

### Reading the initialization files
Once the paths are defined, read the initialization files:
```
config, params, instruments, samples, obs_constraints = io.read_init_files(paths_to_files)
```
There are 5 files in the `init_dir`:
```
init
├── config.yml
├── instruments.yml
├── obs_constraints.yml
├── parameters.yml
└── samples.yml
```
Of these 5 files, you will probably want to modify the `config.yml` and the `parameters.yml` to fit your needs.

### Create the configuration
Once you have read the initialization files, create the configuration for the GRB population with:
```
incl_samples, incl_instruments, incl_constraints = io.create_config(config,
                                                                    samples,
                                                                    instruments,
                                                                    obs_constraints)
```
This function sets up the samples, the instruments and the observational constraints that should be included.

### Setting up the cosmology
The last step before creating the GRB population is to set up the cosmology.
If you want to use the default cosmology do:
```
from cosmology import init_cosmology
cosmo = init_cosmology(paths_to_dir['cosmo'])
```
If not, to change the cosmology, simply do:
```
from cosmology import create_cosmology
cosmo = create_cosmology(OmegaM=0.3, OmegaL=0.7, h=0.7, z_step=0.001, zmax=100)
```

### Generating the GRB population
Now that all of the initialization stuff is out of the way, create the GRB population with:
```
gp = create_GRB_population_from(Nb_GRBs=config['Nb_GRBs'],
                                cosmo=cosmo,
                                params=params,
                                incl_samples=incl_samples,
                                incl_instruments=incl_instruments,
                                incl_constraints=incl_constraints,
                                output_dir=paths_to_dir['output'],
                                run_mode='debug',
                                savefig=True)
```

