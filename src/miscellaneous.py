import logging
import yaml
import numpy as np
import pandas as pd
import argparse

log = logging.getLogger(__name__)


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


def included_instruments(incl_samples, instruments):
    """
        Create the list of the instruments which will be needed by the
        code, given the samples required by the user.
    """
    incl_instruments = {incl_samples[s]['instrument']:
                        instruments[incl_samples[s]['instrument']] for s in incl_samples}

    log.info(f"Including instruments: {list(incl_instruments.keys())}")
    log.debug("Including instruments:\n" + str(yaml.dump(incl_instruments, indent=4)))

    return incl_instruments


def log_to_lin(log_x, log_x_errp, log_x_errm=None):
    """
        Takes logscale data with errors and converts to linear scale with correct error propagation.
        If log_x_errm is not provided, errors are assumed symmetric.
        Returns : x, x_errp, x_errm
    """
    if log_x_errm is None:
        log_x_errm = log_x_errp
    x = 10**log_x
    x_errp = x * (10**log_x_errp - 1.0)
    x_errm = x * (1.0 - 10**(-log_x_errm))

    return x, x_errp, x_errm


def lin_to_log(x, x_errp, x_errm=None):
    """
        Takes linear data with errors and converts to logscale with correct error propagation.
        If x_errm is not provided, errors are assumed symmetric.
        Returns : log_x, log_x_errp, log_x_errm
    """
    if x_errm is None:
        x_errm = x_errp
    log_x = np.log10(x)
    log_x_errp = np.log10((x + x_errp) / x)
    log_x_errm = np.log10(x / (x - x_errm))

    return log_x, log_x_errp, log_x_errm


def lin_to_log_ndarray(x):
    """
    Helper function that turns an ndarray in linscale to logscale with proper error propagation.
    Assumes the data is in the form of the output of read_column (ndarray).
    i.e. :  x[0] = data (float)
            x[1] = error plus (float)
            x[2] = error minus (float)
            x[3] = upper limit (bool)
            x[4] = lower limit (bool)
    """

    log_x = np.zeros(x.shape)
    log_x[0] = np.log10(x[0])
    log_x[1] = np.log10(x[0] + x[1]) - log_x[0]
    log_x[2] = log_x[0] - np.log10(x[0] - x[2])
    log_x[3] = x[3]
    log_x[4] = x[4]

    return log_x


def filter_df(df, filtering_key, lim_min=None, lim_max=None, equal=None, errors='raise', strip=None, string=False):
    """
        Filter a df using a criteria based on filtering_key
    """
    if not string:
        df[filtering_key] = pd.to_numeric(df[filtering_key], errors=errors)
    else:
        df[filtering_key] = df[filtering_key].str.strip(strip)

    if (lim_min is None) and (lim_max is None) and (equal is None):
        raise ValueError
    elif (lim_min is not None) and (lim_max is not None) and (equal is None):
        cond_min = df[filtering_key] >= lim_min
        cond_max = df[filtering_key] <= lim_max
        cond = cond_min & cond_max
    elif lim_min is not None:
        cond = df[filtering_key] >= lim_min
    elif lim_max is not None:
        cond = df[filtering_key] <= lim_max
    elif equal is not None:
        cond = df[filtering_key] == equal
    else:
        print("lim_min = {}".format(lim_min))
        print("lim_max = {}".format(lim_max))
        print("equal = {}".format(equal))
        raise ValueError("You cannot provide all these arguments. Filtering must be "
                         "lim_min and/or lim_max OR equal")
    df_out = df[cond].copy()
    return df_out


def mask_ndarray(ndarray, mask):
    """
        Helper function to easily mask a ndarray output from read_data
    """
    if len(mask) != ndarray.shape[1]:
        print("[Error] in mask_ndarray : mask and array length are not the same")
    masked_ndarray = []
    for i in range(ndarray.shape[0]):
        masked_ndarray.append(ndarray[i,mask])
    masked_ndarray = np.asarray(masked_ndarray)

    return masked_ndarray


def chi2_func(mod, hist, err):

    return np.sum((np.abs(mod-hist)/err)**2)


def str2bool(s):
    if s.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif s.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected')
