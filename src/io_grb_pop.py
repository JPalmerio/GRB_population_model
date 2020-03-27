import logging
import yaml
import datetime
import pickle
import numpy as np
import pandas as pd
from pathlib import Path

log = logging.getLogger(__name__)

root_dir = Path(__file__).resolve().parents[1]


def load_GRBPopulation_from(fname):
    """
        Load a population from a pickle file
    """
    with open(fname, 'rb') as f:
        GRB_pop = pickle.load(f)
    log.debug(GRB_pop.summary())
    return GRB_pop


def generate_paths(conf_fname=None, param_fname=None, init_dir=None):
    """
        Generate the various paths for the necessary input and output
        files.
    """
    log.debug(f"Root directory is {root_dir}")

    obs_dir = root_dir/'observational_constraints'
    output_dir = root_dir/'model_outputs'
    data_dir = root_dir/'data'
    cosmo_dir = data_dir/'cosmology'
    ECLAIRs_dir = data_dir/'ECLAIRs'
    if init_dir is None:
        init_dir = root_dir/'init'

    # Input files
    if conf_fname is None:
        config_file = init_dir/'config.yml'
    else:
        config_file = init_dir/conf_fname
    if param_fname is None:
        param_file = init_dir/'parameters.yml'
    else:
        param_file = init_dir/param_fname

    instrum_file = init_dir/'instruments.yml'
    sample_file = init_dir/'samples.yml'
    obs_constraints_file = init_dir/'obs_constraints.yml'

    paths_to_dir = {'root': root_dir,
                    'obs': obs_dir,
                    'init': init_dir,
                    'output': output_dir,
                    'data': data_dir,
                    'cosmo': cosmo_dir,
                    'ECLAIRs': ECLAIRs_dir}

    paths_to_files = {'config': config_file,
                      'param': param_file,
                      'instrum': instrum_file,
                      'sample': sample_file,
                      'obs_constraints': obs_constraints_file}

    paths_to_dir_str = f"""
    'root'    : {root_dir}
    'obs'     : {obs_dir}
    'init'    : {init_dir}
    'output'  : {output_dir}
    'data'    : {data_dir}
    'cosmo'   : {cosmo_dir}
    'ECLAIRs' : {ECLAIRs_dir}"""

    paths_to_files_str = f"""
    'config'          : {config_file}
    'param'           : {param_file}
    'instrum'         : {instrum_file}
    'sample'          : {sample_file}
    'obs_constraints' : {obs_constraints_file}"""

    log.debug("Directory paths :" + paths_to_dir_str)
    log.debug("File paths :" + paths_to_files_str)

    return paths_to_dir, paths_to_files


def read_init_files(paths_to_files):
    """
        Read the various files that set up the configuration for the
        code.
    """

    with open(paths_to_files['config'], 'r') as f:
        config = yaml.safe_load(f)
    with open(paths_to_files['param'], 'r') as f:
        params = yaml.safe_load(f)
    with open(paths_to_files['instrum'], 'r') as f:
        instruments = yaml.safe_load(f)
    with open(paths_to_files['sample'], 'r') as f:
        samples = yaml.safe_load(f)
    with open(paths_to_files['obs_constraints'], 'r') as f:
        obs_constraints = yaml.safe_load(f)

    log.debug("Input configuration:\n" + str(yaml.dump(config, indent=4)))
    log.debug("Input parameters:\n" + str(yaml.dump(params, indent=4)))

    return config, params, instruments, samples, obs_constraints


def create_config(config, samples, instruments, obs_constraints):
    """
        Filter the samples, instruments, and observation constraint
        dictionaries to return only the ones included in the config.
        The samples, instruments and obs_constraints inputs represent
        all possible choices, and the specific ones to be used are
        specified by the user in the config file.
    """
    incl_samples = included_samples(config['samples'], samples)
    incl_instruments = included_instruments(incl_samples, instruments)
    incl_constraints = included_constraints(config['constraints'], obs_constraints)
    incl_constraints = load_observational_constraints(incl_constraints)
    return incl_samples, incl_instruments, incl_constraints


def included_samples(config_samples, samples):
    """
        Create the list of the sample names which will be used by the
        code.
    """
    try:
        incl_samples = {sample_name: samples[sample_name] for sample_name in config_samples}
    except KeyError as e:
        raise ValueError(f"Sample {e} does not exist. Possible samples are {list(samples.keys())}")

    log.info(f"Including samples: {list(incl_samples.keys())}")
    log.debug("Including samples:\n" + str(yaml.dump(incl_samples, indent=4)))

    return incl_samples


def included_constraints(config_constr, constraints):
    """
        Create the list of the instruments which will be needed by the
        code, given the samples required by the user.
    """
    try:
        incl_constraints = {constr_name: constraints[constr_name] for constr_name in config_constr}
    except KeyError as e:
        raise ValueError(f"Constraint {e} does not exist."
                         f"Possible choices are {list(constraints.keys())}")

    log.info(f"Including constraints: {list(incl_constraints.keys())}")
    log.debug("Including constraints:\n" + str(yaml.dump(incl_constraints, indent=4)))

    return incl_constraints


def included_instruments(incl_samples, instruments):
    """
        Create the list of the instruments which will be needed by the
        code, given the samples required by the user.
    """
    incl_instruments = {incl_samples[s]['instrument']: instruments[incl_samples[s]['instrument']]
                        for s in incl_samples}

    log.info(f"Including instruments: {list(incl_instruments.keys())}")
    log.debug("Including instruments:\n" + str(yaml.dump(incl_instruments, indent=4)))

    return incl_instruments


def load_observational_constraints(obs_constraints):
    """
        Read the observational data once and for all and store them in
        a dictionary to be used by the fitting procedure.
    """

    for name, constraint in obs_constraints.items():
        bins, hist, err = read_constraint(name, last_bin_edge=constraint['last_bin'])
        constraint['bins'] = bins
        constraint['hist'] = hist
        constraint['err'] = err
    return obs_constraints


def read_constraint(name, last_bin_edge, bins_log=False, density=False, verbose=False):
    """
        A convenience function to read the observational constraint from
        a file. Assumes the format of the file is:
        1. left bin edge 2. histogram of counts 3. error
        You must provide a last bin edge as well.
    """

    valid_names = ['Stern', 'EpGBM', 'eBAT6', 'Kommers']

    if name not in valid_names:
        raise ValueError('Constraint name must be one of {}'.format(valid_names))

    fname = root_dir/f'observational_constraints/{name}.txt'

    bins = read_column(fname, 0, array=False)
    hist = read_column(fname, 1)
    err = read_column(fname, 2)
    bins.append(last_bin_edge)  # append right edge of last bin
    bins = np.array(bins)

    if bins_log:
        bins = np.log10(bins)

    if density:
        N_EpGBM = hist.sum()
        delta_bin = bins[1:]-bins[:-1]
        hist /= (N_EpGBM * delta_bin)
        err /= (N_EpGBM * delta_bin)

    if verbose:
        ln_oi = 0.
        for i, val in enumerate(hist):
            ln_oi += val*np.log(val) - val
        print(f"ln(o_i!) = {ln_oi} from EpGBM histogram")
    return bins, hist, err


def read_logRlogN(fname=None):
    """
        Create the histogram from [Stern et al. 2001]
        (https://ui.adsabs.harvard.edu/abs/2001ApJ...563...80S/abstract)
        Fig. 23 rebinned to avoid low counts
        This histogram is log( Rate/delta(log pflx) ) in units of
        [GRBs/yr/log pflx in 4 pi] (i.e. corrected for solid angle and
        live time of the search):
        <Omega_BATSE> = 0.7 * 4pi
        T_live_BATSE = 6.25 yr (T_tot_BATSE = 9.1 yr) see Goldstein 2013
    """
    if fname is None:
        fname = root_dir/'resources/Stern_logRlogN.txt'

    bins = read_column(fname, 0, array=False)
    hist = read_column(fname, 1)
    err = read_column(fname, 2)
    bins.append(50.0)  # append right edge of last bin
    bins = np.array(bins)

    return bins, hist, err


def create_output_dir(paths_to_dir, dir_name, overwrite=False):
    """
        Create the output directory where the results from the code will
        be saved. This directory is named as YYMMDD_{output_dir}
        with YYMMDD being the date and {output_dir} defined in the
        file config.yml
    """

    today = datetime.datetime.now().strftime('%y%m%d')
    _output_dir = paths_to_dir['output'] / '_'.join([today, dir_name])
    if overwrite:
        # if in overwrite mode, write over any old outputs
        if not _output_dir.is_dir():
            _output_dir.mkdir()
    else:
        # otherwise, make sure directory doesn't exist, if it does
        # append _n at the end and create a new one to avoid overwriting
        n = 2
        while _output_dir.is_dir():
            log.warning(f"{_output_dir} already exists...")
            _output_dir = paths_to_dir['output'] / '_'.join([today, dir_name, f'{n}'])
            log.info(f"Trying {_output_dir}")
            n += 1
            if n > 10:
                raise ValueError('Over 10 directories already exist. Something probably went wrong')
        _output_dir.mkdir()
    paths_to_dir['output'] = _output_dir

    log.info(f"Output directory updated to {_output_dir}")

    return


def read_column(filename, column_nb, end=None, dtype=float, array=True, splitter=None, stripper=None, verbose=False):
    """
    Function used to read ASCII files.
    It will skip lines starting with '#', '!' or '%'.

    Parameters
    ----------
    filename : [str]
        Name of the file containing the data

    column_nb: [int]
        Number of the column for the data.

    end   : [int]
        Number of lines to read. Note: commented lines (i.e. starting with '#', '!', or '%') do not count as lines for this purpose.
        Default is None, which reads the whole file.

    dtype : [data-type]
        Type of the returned data. Default is float.

    array : [boolean]
        If True returns xdata as an array rather than a list. Default is True (arrays are faster)

    splitter : [str]
        String to use as a delimiter between columns. Default is None (uses default for str.split() which is a whitespace)

    stripper : [str]
        String to strip at the end of each line. Default is None (uses default for str.strip() which is a whitespace)

    Returns
    -------
    xdata : [array/list]
        x data
    """

    nan = False

    xdata = []
    i_counter = 0
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if len(line) != 0:
                if line[0] != '#' and line[0] != '!' and line[0] != '%':
                    i_counter += 1
                    line = line.strip(stripper)
                    columns = line.split(splitter)
                    try:
                        xdata.append(dtype(columns[column_nb]))
                    except ValueError:
                        nan = True
                        xdata.append(np.nan)
                        if verbose:
                            log.error("In read_column for {}: could not convert {} to {}, \
                                       so added NaN.".format(filename, columns[column_nb], dtype))
                    except IndexError:
                        if verbose:
                            log.error("In read_column for {}: no data found for column {:d}, line {:d}. \
                                       Input will be NaN.".format(filename, column_nb, i))
            if (end is not None) and (i_counter >= end):
                break
    if array:
        xdata = np.asarray(xdata, dtype=dtype)
    if nan:
        log.warning("In read_column for {}: some strings could not be converted to {} \
            NaNs were added instead.".format(filename, dtype))

    return xdata


def read_data(filename, column_nb, end=None, err=True, dtype=float, single_err=False, splitter=None, stripper=None, verbose=False):
    """
    Helper function that reads a table and returns the data, with errors and upper or lower limits.
    If err is True it assumes format of file is :     1) data     2) error plus     3) error minus
    If err is False, it sets them to 0.
    If single_err is True it assumes format of file is : 1) data     2) error
    If there is no data it fills the errors with NaN.
    Otherwise if the data exists but it doesn't find errors it will set them to 0.
    Upper and lower limits and converted to 1 for True, and 0 for False in the requested dtype.
    Returns data as a list of numpy.ndarray with by default:
        data[0] = data        (float)
        data[1] = plus error  (float)
        data[2] = minus error (float)
        data[3] = upper limit (float)
        data[4] = lower limit (float)
    See help for read_column for more information on arguments.
    """
    data = [[], [], [], [], []]
    f = open(filename, 'r')
    i = 1
    i_counter = 0
    if err is True:
        for line in f:
            if line[0] != '#' and line[0] != '!' and line[0] != '%':
                i_counter += 1
                line = line.strip(stripper)
                columns = line.split(splitter)
                try:
                    # Check for upper limit
                    if '<' in columns[column_nb]:
                        # Remove '<' from data
                        place = columns[column_nb].find('<')
                        columns[column_nb] = columns[column_nb][place+1:]
                        # Append upper limit
                        data[3].append(True)
                        # Set errors to 0.0, and lower limit to false
                        data[1].append(0.)  # error plus
                        data[2].append(0.)  # error minus
                        data[4].append(False)

                    # Check for lower
                    elif '>' in columns[column_nb]:
                        # Remove '<' from data
                        place = columns[column_nb].find('>')
                        columns[column_nb] = columns[column_nb][place+1:]
                        # Append lower limit
                        data[4].append(True)
                        # Set errors to 0.0, and upper limit to false
                        data[1].append(0.)  # error plus
                        data[2].append(0.)  # error minus
                        data[3].append(False)

                    # if no lower or upper limits, append errors
                    else:
                        data[3].append(False)
                        data[4].append(False)
                        try:
                            data[1].append(float(columns[column_nb+1]))  # error plus
                        except IndexError:
                            data[1].append(None)  # error plus
                        except ValueError:
                            data[1].append(0.)  # error plus
                        try:
                            data[2].append(np.abs(float(columns[column_nb+2])))  # error minus
                        except IndexError:
                            data[2].append(None)  # error minus
                        except ValueError:
                            data[2].append(0.)  # error minus

                    # Append data
                    try:
                        data[0].append(dtype(columns[column_nb]))
                    except ValueError:
                        if verbose:
                            print("[Warning] in read_data : Couldn't convert %s to %s for column %d, line %d in file %s. Replacing with NaN." % (columns[column_nb], dtype,column_nb, i, filename))
                        data[0].append(np.nan)
                # If no data
                except IndexError:
                    if verbose:
                        print('[Warning] in read_data : No data found for column %d, line %d in file %s. Input will be NaN.' % (column_nb, i, filename))
                    data[0].append(np.nan)  # data
                    data[1].append(np.nan)  # error plus
                    data[2].append(np.nan)  # error minus
                    data[3].append(np.nan)  # upper limit
                    data[4].append(np.nan)  # lower limit
            i += 1
            if (end is not None) and (i_counter >= end):
                break
    else:
        for line in f:
            if line[0] != '#' and line[0] != '!' and line[0] != '%':
                i_counter += 1
                line = line.strip(stripper)
                columns = line.split(splitter)
                try:
                    # Check for upper limit
                    if '<' in columns[column_nb]:
                        # Remove '<' from data
                        place = columns[column_nb].find('<')
                        columns[column_nb] = columns[column_nb][place+1:]
                        # Append upper limit
                        data[3].append(True)
                        # Set lower limit to false
                        data[4].append(False)

                    # Check for lower
                    elif '>' in columns[column_nb]:
                        # Remove '<' from data
                        place = columns[column_nb].find('>')
                        columns[column_nb] = columns[column_nb][place+1:]
                        # Append lower limit
                        data[4].append(True)
                        # Set upper limit to false
                        data[3].append(False)

                    else:
                        data[3].append(False)
                        data[4].append(False)

                    # Append data
                    try:
                        data[0].append(dtype(columns[column_nb]))
                    except ValueError:
                        if verbose:
                            print("[Warning] in read_data : Couldn't convert %s to %s for column %d, line %d in file %s. Replacing with NaN." % (columns[column_nb], dtype,column_nb, i, filename))
                        data[0].append(np.nan)
                    data[1].append(0.0)  # error plus
                    data[2].append(0.0)  # error minus

                except IndexError:
                    # print 'Warning : No data found for column %d, line %d in file %s. Input will be None.' %(column_nb, i, filename)
                    data[0].append(np.nan)  # data
                    data[1].append(np.nan)  # error plus
                    data[2].append(np.nan)  # error minus
                    data[3].append(np.nan)  # upper limit
                    data[4].append(np.nan)  # lower limit
            if (end is not None) and (i_counter >= end):
                break
    if single_err:
        data[2] = data[1]

    data[0] = np.array(data[0]).astype(dtype)
    data[1] = np.array(data[1]).astype(dtype)
    data[2] = np.array(data[2]).astype(dtype)
    data[3] = np.asarray(data[3]).astype(dtype)
    data[4] = np.asarray(data[4]).astype(dtype)
    data = np.asarray(data)
    return data


def read_SHOALS_file(SHOALS_file=None):
    if SHOALS_file is None:
        SHOALS_file = root_dir/'catalogs/SHOALS_cat/SHOALS_cat.txt'

    df_obs = pd.read_csv(SHOALS_file, sep='\t', header=3, low_memory=False)
    keys = ['S_BAT', 'z']
    for key in keys:
        df_obs[key] = pd.to_numeric(df_obs[key], errors='coerce')
    df_obs = df_obs.dropna()
    return df_obs
