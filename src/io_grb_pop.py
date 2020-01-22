import logging
import yaml
import datetime
import numpy as np
from pathlib import Path

log = logging.getLogger(__name__)

root_dir = Path(__file__).resolve().parents[1]


def generate_paths():
    """
        Generate the various paths for the necessary input and output
        files.
    """
    log.debug(f"Root directory is {root_dir}")

    obs_dir = root_dir/'observational_constraints'
    init_dir = root_dir/'init'
    output_dir = root_dir/'model_outputs'
    data_dir = root_dir/'data'
    cosmo_dir = data_dir/'cosmology'
    ECLAIRs_dir = data_dir/'ECLAIRs'

    # Input files
    configfile = init_dir/'config.yml'
    paramfile = init_dir/'parameters.yml'
    instrumfile = init_dir/'instruments.yml'
    samplefile = init_dir/'samples.yml'

    paths_to_dir = {'root': root_dir,
                    'obs': obs_dir,
                    'init': init_dir,
                    'output': output_dir,
                    'data': data_dir,
                    'cosmo': cosmo_dir,
                    'ECLAIRs': ECLAIRs_dir}

    paths_to_files = {'config': configfile,
                      'param': paramfile,
                      'instrum': instrumfile,
                      'sample': samplefile}

    paths_to_dir_str = f"""
    'root'    : {root_dir}
    'obs'     : {obs_dir}
    'init'    : {init_dir}
    'output'  : {output_dir}
    'data'    : {data_dir}
    'cosmo'   : {cosmo_dir}
    'ECLAIRs' : {ECLAIRs_dir}"""

    paths_to_files_str = f"""
    'config'  : {configfile}
    'param'   : {paramfile}
    'instrum' : {instrumfile}
    'sample'  : {samplefile}"""

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

    log.debug("Input configuration:\n" + str(yaml.dump(config, indent=4)))
    log.debug("Input parameters:\n" + str(yaml.dump(params, indent=4)))

    return config, params, instruments, samples


def create_output_dir(paths_to_dir, dir_name, run_mode=None):
    """
        Create the output directory where the results from the code will
        be saved. This directory is named as YYMMDD_{output_dir}
        with YYMMDD being the date and {output_dir} defined in the
        file config.yml
    """

    today = datetime.datetime.now().strftime('%y%m%d')
    _output_dir = paths_to_dir['output'] / '_'.join([today, dir_name])
    if run_mode == 'debug':
        if not _output_dir.is_dir():
            _output_dir.mkdir()
    else:
        n = 2
        while _output_dir.is_dir():
            log.warning(f"{_output_dir} already exists...")
            _output_dir = paths_to_dir['output'] / '_'.join([today, dir_name, f'{n}'])
            log.info(f"Trying {_output_dir}")
            n += 1
            if n > 100:
                raise ValueError('Over 100 directories already exist. Something probably went wrong')
        _output_dir.mkdir()
    paths_to_dir['output'] = _output_dir

    log.info(f"Output directory updated to {_output_dir}")

    return


def read_column(filename, column_nb, end=None, dtype=float, array=True,
                splitter=None, stripper=None, verbose=False):
    """
    Function used to read ASCII (or fits) files.
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
    i = 0
    f = open(filename, 'r')
    for line in f:
        if len(line) != 0:
            if line[0] != '#' and line[0] != '!' and line[0] != '%':
                i_counter += 1
                if stripper is not None:
                    line = line.strip(stripper)
                else:
                    line = line.strip()
                if splitter is not None:
                    columns = line.split(splitter)
                else:
                    columns = line.split()
                try:
                    xdata.append(dtype(columns[column_nb]))
                except ValueError:
                    nan = True
                    xdata.append(np.nan)
                except IndexError:
                    if verbose:
                        log.warning('In read_data : No data found for column %d, line %d in file %s. Input will be NaN.' % (column_nb, i, filename))
        i += 1
        if (end is not None) and (i_counter >= end):
            break
    if(array):
        xdata = np.asarray(xdata, dtype=dtype)
    f.close()
    if nan:
        if verbose:
            log.warning("In read_column for %s. Could not convert string to %s, so added NaN." % (filename, dtype))

    return xdata
