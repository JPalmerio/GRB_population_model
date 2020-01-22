import logging
import yaml
import numpy as np
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


def chi2_func(mod, hist, err):

    return np.sum((np.abs(mod-hist)/err)**2)


def str2bool(s):
    if s.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif s.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected')
