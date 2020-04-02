# GRB Population model

`grbpop` allows one to create a population of Gamma-Ray Bursts (GRBs).
Using Monte Carlo sampling from input distributions, one can then compare the predictions for various observational constraints.
If you have any questions, please contact <palmerio@iap.fr>
# Installation
There is no pip installation for now, the simplest way to get this code is to git clone this repository:
```
git clone https://github.com/JPalmerio/GRB_population_model.git
```

# Basic usage
The arborescence is important for this software.
Make sure you have the following directory structure (this is the default if you clone the directory):
```bash
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
The `data/cosmology` directory is also optional; you can create your own cosmology if you wish (see ***Setting up the cosmology***).

A basic example is shown below, if you want more details see ***Detailed usage*** below:
```python
import sys
src_dir = 'grbpop/'
sys.path.insert(0, src_dir)
from cosmology import init_cosmology
from GRB_population import create_GRB_population_from
import io_grb_pop as io
import numpy as np
import logging

log = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(asctime)s.%(msecs)03d [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
logging.getLogger('matplotlib').setLevel(logging.WARNING)


# Define the paths used by the code
paths_to_dir, paths_to_files = io.generate_paths(conf_fname='config_simple_example.yml',
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
```

That's it. Everything you need is in the `gp` object which is an instance of the `GRBPopulation` class.



# Detailed usage
You will want to start your code with the following:
```python
import io_grb_pop as io
from GRB_population import GRBPopulation, create_GRB_population_from
import logging
log = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG,
                    format='%(asctime)s.%(msecs)03d [%(levelname)s] %(message)s',
                    datefmt='%H:%M:%S')
```
You can change the logging level to `INFO` if you wish.

## Setting up the paths

The first thing to do is set up the paths to the various files so that the code knows where to look for what:
```python
paths_to_dir, paths_to_files = io.generate_paths(init_dir=None)
```
If you did not follow the directory structure presented above, you can provide a `init_dir` argument to tell the code where to look for the initialization files.

## Reading the initialization files

Once the paths are defined, read the initialization files:
```python
config, params, instruments, samples, obs_constraints = io.read_init_files(paths_to_files, obs_dir=None)
```
There are 5 files in the `init_dir`:
```bash
.
├── init
│   ├── config.yml
│   ├── parameters.yml
│   ├── instruments.yml
│   ├── samples.yml
│   └── obs_constraints.yml
└── observational_constraints
    ├── EpGBM.txt
    ├── Stern.txt
    └── eBAT6.txt
```
Of these 5 files, you will probably want to modify the `config.yml` and the `parameters.yml` to fit your needs.
You can also specify where the observational constraints are to be found by passing an `obs_dir` argument; by default the above structure is assumed.
If you want to add observational constraints see [*this*]().

## Create the configuration

Once you have read the initialization files, create the configuration for the GRB population with:
```python
incl_samples, incl_instruments, incl_constraints = io.create_config(config,
                                                                    samples,
                                                                    instruments,
                                                                    obs_constraints)
```
This function sets up the samples, the instruments and the observational constraints that should be included.

## Setting up the cosmology

The last step before creating the GRB population is to set up the cosmology.
If you want to use the default cosmology do:
```python
from cosmology import init_cosmology
cosmo = init_cosmology(paths_to_dir['cosmo'])
```
If not, to change the cosmology, simply do:
```python
from cosmology import create_cosmology
cosmo = create_cosmology(OmegaM=0.3, OmegaL=0.7, h=0.7, z_step=0.001, zmax=100)
```

## Generating the GRB population

Now that all of the initialization stuff is out of the way, creating the GRB population is as simple as:
```python
io.create_output_dir(paths_to_dir=paths_to_dir, dir_name=config['output_dir'], overwrite=True)

# Create the class instance
gp = GRBPopulation(Nb_GRBs=config['Nb_GRBs'], output_dir=paths_to_dir['output'])

# Draw the properties : z, L, Ep, alpha, beta, t90, Cvar
gp.draw_GRB_properties(cosmo=cosmo, params=params, run_mode='debug', savefig=True)

# Calculate the peak photon and energy fluxes for the included instruments
# Also calculate the probability of detection for the included samples
gp.calculate_quantities(instruments=incl_instruments, samples=incl_samples)

# Create the mock constraints for the included constraints
gp.create_mock_constraints(constraints=incl_constraints)

# Compare the mock and real constraints and calculate the likelihood of the population
gp.compare_to_observational_constraints(constraints=incl_constraints)

# Normalize the population to the Intensity constraint of Stern et al. 2001
gp.normalize_to_Stern()

print(gp.summary())
```
This might take a while depending on `Nb_GRBs`; I recommend $10^5$ for fast but reliable numbers.
From this point on, you can access all the aspects of the population through the `gp` object.

### Properties of the GRB population

To get the properties of each GRB, simply do:
```python
df = gp.properties
```
This returns a [pandas](https://pandas.pydata.org/docs/getting_started/10min.html) `DataFrame` of length `Nb_GRBs` and with columns depending on your desired configuration. Here is an example:
```
          z           D_L             L            Ep     alpha     beta  \
0      1.432  10391.248459  1.532278e+50    269.700917  0.883468  2.75087   
1      3.170  27630.627622  3.344103e+50    211.470520  1.057890  2.49724   
2      1.791  13726.341548  8.953720e+51   1613.322041  0.509119  2.19460   
3      2.447  20166.409982  5.994049e+50    352.916938 -0.227648  2.56481   
4      0.564   3267.229951  5.779901e+50    476.598862  0.855250  2.41085   
...      ...           ...           ...           ...       ...      ...   
99995  1.715  13007.160092  6.000251e+52  14978.008833  0.631182  2.88288   
99996  1.355   9698.523642  2.735302e+52   2674.374924  0.328647  2.68191   
99997  1.675  12631.304370  4.605083e+50    150.124012 -0.054855  2.14142   
99998  4.192  38634.269398  6.952059e+50    604.285453  0.915621  3.14171   
99999  2.150  17203.801554  3.449701e+50    138.688225  0.707670  4.04809   

          ktild        Epobs  pht_pflx_BATSE        t90      t90obs      Cvar  \
0      0.954909   110.896759        0.034102  84.331929  205.095252  0.241399   
1      0.637547    50.712355        0.007777   6.495532   27.086367  0.214844   
2      0.680962   578.044443        0.171092  19.495005   54.410560  0.315690   
3      3.184627   102.383794        0.036779   8.250524   28.439558  0.144807   
4      0.754715   304.730730        0.661843  12.631445   19.755580  0.131637   
...         ...          ...             ...        ...         ...       ...   
99995  1.409505  5516.762001        0.203914   3.885583   10.549357  0.500129   
99996  1.852285  1135.615679        0.834032   3.111309    7.327134  0.395769   
99997  0.953798    56.121126        0.028595  30.114104   80.555228  0.049988   
99998  1.025939   116.387799        0.012345   7.425105   38.551145  0.326739   
99999  1.503101    44.028008        0.025734   3.933657   12.391019  0.625070   

               Eiso  erg_pflx_BATSE  pht_flnc_BATSE  erg_flnc_BATSE  \
0      3.119362e+51    5.747455e-09        1.688369    2.845557e-07   
1      4.666793e+50    1.199519e-09        0.045259    6.980426e-09   
2      5.510463e+52    3.948043e-08        2.938827    6.781508e-07   
3      7.161282e+50    6.222786e-09        0.151465    2.562693e-08   
4      9.610598e+50    1.347467e-07        1.721161    3.504167e-07   
...             ...             ...             ...             ...   
99995  1.166024e+53    4.948661e-08        1.075859    2.610932e-07   
99996  3.368141e+52    2.083402e-07        2.418571    6.041558e-07   
99997  6.932280e+50    4.758908e-09        0.115146    1.916327e-08   
99998  1.686621e+51    2.091221e-09        0.155502    2.634138e-08   
99999  8.482159e+50    3.288671e-09        0.199320    2.547159e-08   

       pdet_Stern  
0             0.0  
1             0.0  
2             1.0  
3             0.0  
4             1.0  
...           ...  
99995         1.0  
99996         1.0  
99997         0.0  
99998         0.0  
99999         0.0  

[100000 rows x 17 columns]
```

#### Names of the columns and units

The column names for the `properties` follow a certain structure. There are the properties that are drawn for each LGRB:

| Key | Units | Description |
| :--- | :---: | :--- |
| z | - | The redshift of the LGRB source
| L | erg/s | The isotropic-equivalent bolometric peak luminosity
| Ep | keV | The peak energy of the photon spectrum in the source frame
| alpha | - | The low-energy slope of the photon spectrum $L_E/E$ in the source frame
| beta | - | The high-energy slope of the photon spectrum $L_E/E$ in the source frame
| t90 | s | The duration over which 5 to 95 \% of photons are emitted in the source frame
| Cvar | - | The variability coefficient defined as the mean luminosity divided by the peak luminosity

And the properties that are calculated:
| Key | Units | Description |
| :--- | :---: | :--- |
| D_L | Mpc | The luminosity distance to the LGRB source
| Eiso | erg | The isotropic-equivalent bolometric energy
| ktild | - | Normalization factor for the photon spectrum, given alpha and beta
| Epobs | keV | The peak energy of the photon spectrum in the observer frame
| t90obs | s | The duration over which 5 to 95 \% of photons are received in the observer frame
| pht_pflx_*instr* | ph/s/cm2 |  The peak photon flux for a given *instrument*
| erg_pflx_*instr* | erg/s/cm2 |  The peak energy flux for a given *instrument*
| pht_flnc_*instr* | ph/cm2 |  The photon fluence for a given *instrument*
| erg_flnc_*instr* | erg/cm2 |  The energy fluence for a given *instrument*
| pdet_*sample* | - | The probability of detection for a given *sample*

##### Special case of *ECLAIRs*

For *ECLAIRs*, the detection is more complicated as there is a peak flux mode and a fluence mode.
For this reason there is a `pdet_ECLAIRs_tot`, for the total probability of detection, a `pdet_ECLAIRs_pht_cts` for the photon count probability of detection and a `pdet_ECLAIRs_pht_flnc` for a photon fluence probability of detection.
In the same way, the `pflx` prefix is replaced by `cts` in the table above.

#### Samples and instruments

Each sample[^1] is defined by a minimum peak photon flux and an instrument.
[^1]: Except the ECLAIRs and SHOALS samples but they had to be hard coded.
You can create a new sample by simpling modifying the `samples.yml` file in the `init` directory:
```yml
YourSampleName:
  instrument: 'YourInstrument'
  pflx_min: 0.42     # ph/s/cm2 in YourInstrument's energy band
```
For the moment the code only supports samples defined by a peak photon flux cut.
If your instrument is not among the instruments already defined in the `instruments.yml` file, you can add it as:
```yml
YouInstrument:
  Emin: 10.  # keV
  Emax: 1000. # kev
```

To include this sample in the GRB population simply add it to the list of samples in the `config.yml` file:
```yml
samples: ['Stern', 'EpGBM', 'eBAT6', 'YourSampleName']
```

#### Observational constraints

The observational constraints are assumed to be histograms of number counts for a given instrument.
The data is collected in the directory `observational_constraints` in the form of a text file with 3 columns (`left bin edge`, `histogram`, `error`) separated by tabs.
The name of the file should be the name of the constraint as defined in the `init/obs_constraints.yml` file.
To add a new constraint, add to the `obs_constraints.yml`:
```yml
YourConstraint:
  instrument: 'YourInstrument'
  val_min: 0.42     # ph/s/cm2 in YourInstrument's energy band
  prop_min: 'pht_pflx_YourInstrument'
  quantity: 'z'
  sum_ln_oi_factorial: 6870.371682542675 # this is the sum over all bins of ln(o_i!) where o_i in the number count in bin i
  last_bin: 6.0
```

You can access the mock constraints through:
```python
gp.mock_constraints
```
This returns a dictionary with each key being the name of the constraint as defined in the `config.yml`:
```yml
constraints: ['Stern', 'EpGBM', 'eBAT6', 'YourConstraint']
```
You can access the content of the dictionary with:
```python
>>> gp.mock_constraints['Stern']

{'val_min': 0.066825,
 'prop_min': 'pht_pflx_BATSE',
 'quantity': 'pht_pflx_BATSE',
 'bins': array([ 0.066825 ,  0.0841276,  0.10591  ,  0.133333 ,  0.167857 ,
         0.211319 ,  0.266035 ,  0.334918 ,  0.421637 ,  0.53081  ,
         0.66825  ,  0.841276 ,  1.0591   ,  1.33333  ,  1.67857  ,
         2.11319  ,  2.66035  ,  3.34918  ,  4.21637  ,  5.3081   ,
         6.6825   ,  8.41277  , 10.591    , 13.3333   , 16.       ,
        20.       , 28.       , 50.       ]),
 'hist_unnormed': array([3544, 3280, 2903, 2805, 2463, 2133, 2043, 1812, 1591, 1391, 1239,
        1100,  933,  803,  732,  598,  488,  418,  320,  249,  214,  145,
          99,   54,   54,   69,   66]),
 'err_unnormed': array([59.53150426, 57.27128425, 53.87949517, 52.96225071, 49.62862077,
        46.18441296, 45.19955752, 42.56759331, 39.88734135, 37.2961124 ,
        35.19943181, 33.1662479 , 30.5450487 , 28.33725463, 27.05549852,
        24.45403852, 22.09072203, 20.4450483 , 17.88854382, 15.77973384,
        14.62873884, 12.04159458,  9.94987437,  7.34846923,  7.34846923,
         8.30662386,  8.1240384 ]),
 'norm': 0.22923047296917426,
 'hist': array([812.3927962 , 751.87595134, 665.45606303, 642.99147668,
        564.59465492, 488.94859884, 468.31785628, 415.36561702,
        364.70568249, 318.8595879 , 284.01655601, 252.15352027,
        213.87203128, 184.07206979, 167.79670621, 137.07982284,
        111.86447081,  95.8183377 ,  73.35375135,  57.07838777,
         49.05532122,  33.23841858,  22.69381682,  12.37844554,
         12.37844554,  15.81690263,  15.12921122]),
 'err': array([13.64643488, 13.12832358, 12.35082216, 12.14056178, 11.37639221,
        10.58687483, 10.36111595,  9.75778955,  9.14339412,  8.54940548,
         8.0687824 ,  7.60271469,  7.00185596,  6.49576228,  6.20194472,
         5.60561082,  5.06386666,  4.68662809,  4.10059936,  3.61719585,
         3.35335272,  2.76030042,  2.28081441,  1.68449308,  1.68449308,
         1.90413132,  1.86227717])}
```

#### Likelihood

You can see how well your population fit the constraints by accessing the `likelihood_params` attribute:
```python
>>> gp.likelihood_params
{'chi2_Stern': 32.130352211339904,
 'lnL_Stern': -16.065176105669952,
 'lnL_tot': -16.065176105669952,
 'chi2_tot': 32.130352211339904}
```

#### Population normalization
You can find the normalization for your population in:
```python
>>> gp.normalization
{'pseudo_collapse_rate': 2747967987338.383,    # Mpc3
 'T_sim': 28.41981981981982,                   # yr
 'T_sim_err': 1.0241376511646783,              # yr
 'R_intr': 3518.6711468965955,                 # GRB/yr
 'R_intr_err': 126.79896024852597,             # GRB/yr
 'nGRB0': 1.2804629322864484e-09,              # GRB/yr/Mpc
 'nGRB0_err': 4.6142808370682826e-11,          # GRB/yr/Mpc
 }
```
Alternatively, you can provide your own normalization by replacing the following lines:
```python
# Normalize the population to the Intensity constraint of Stern et al. 2001
# gp.normalize_to_Stern()
# Replace with:
gp.normalize_from(nGRB0=1.3e-9, nGRB0_err=0.4e-9) # Units: GRB/yr/Mpc3
```
