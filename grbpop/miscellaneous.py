import logging
import numpy as np
import pandas as pd
from io_grb_pop import root_dir, read_logRlogN

log = logging.getLogger(__name__)


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


def filter_df(df, filtering_key, lim_min=None, lim_max=None, equal=None,
    errors='raise', strip=None, string=False):
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


def calc_cat_duration(fname, verbose=False):
    from datetime import datetime
    data = pd.read_csv(fname, sep='|', header=2, low_memory=False)
    data.rename(columns=lambda x:x.strip(), inplace=True)
    first_date = datetime.strptime(data['name'].min().strip('GRB')[:-3], '%y%m%d')
    last_date = datetime.strptime(data['name'].max().strip('GRB')[:-3], '%y%m%d')
    delta_t = (last_date - first_date).days/364.25
    if verbose:
        print('First GRB detected on {}'.format(first_date))
        print('Last GRB detected on {}'.format(last_date))
        print('Duration: {:.3f} years'.format(delta_t))
    return delta_t


def create_filtered_sample(fname, keys, func=None, func_args={}, ax=None, log=False, kde=True,
    header=2, verbose=False, debug=False, errors='raise', **kwargs):
    """
        Convenience function to quickly plot an observed sample from a given file name.
        A function to filter or cut the sample can be passed as func.
    """
    # Read the entire file
    df_obs = pd.read_csv(fname, sep='|', header=header, low_memory=False)
    # Strip the colum names to remove whitespaces
    df_obs.rename(columns=lambda x:x.strip(), inplace=True)
    # Activate debug to check the column names
    if debug:
        for i,col in enumerate(df_obs.columns):
            print(i,col)

    # Apply function to the data
    if func is None:
        for key in keys:
            df_obs[key] = pd.to_numeric(df_obs[key], errors=errors)
    # If func is a list, iterate through the list and apply each function
    elif isinstance(func, list):
        if not isinstance(func_args, list):
            raise ValueError
        for i, func_i in enumerate(func):
            df_obs = func_i(df_obs, **func_args[i])
        for key in keys:
            df_obs[key] = pd.to_numeric(df_obs[key], errors=errors)
    else:
        df_obs = func(df_obs.copy(), **func_args)
        df_obs[key] = pd.to_numeric(df_obs[key], errors=errors)

    if log:
        for key in keys:
            df_obs[key] = pd.to_numeric(df_obs[key], errors=errors)
            df_obs[key] = np.log10(df_obs[key])
    if verbose:
        print("Sample size :{}".format(len(df_obs.dropna())))

    return df_obs


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


def sort_ndarray(ndarray, sorter=None):
    """
        Sorts a ndarray with finite values !! (will fail if there are NaNs)
        If sorter is none, sorts the array with its values.
        Otherwise, uses sorter as an argsort output.
    """
    if sorter is None:
        arr_ind = ndarray[0].argsort()
    else:
        arr_ind = sorter
    sorted_data = np.zeros(ndarray.shape)
    sorted_data[:,:] = ndarray[:,arr_ind]

    return sorted_data


def calc_rel_errors_GBM_band(df_obs):
    # Calculate relative errors on Band parameters for later filtering
    keys = ['pflx_band_epeak_pos_err','pflx_band_epeak', 'pflx_band_ampl_pos_err', 'pflx_band_ampl']
    for k in keys:
        df_obs[k] = pd.to_numeric(df_obs[k], errors='coerce')
    df_obs['pflx_band_epeak_rel_err'] = df_obs['pflx_band_epeak_pos_err']/df_obs['pflx_band_epeak']
    df_obs['pflx_band_ampl_rel_err'] = df_obs['pflx_band_ampl_pos_err']/df_obs['pflx_band_ampl']
    return df_obs


def global_GRB_rate_Stern():
    """
        Calculate the global LGRB rate by correcting the Stern histogram
        with the efficiency correction [Stern et al. 2001]
        (https://ui.adsabs.harvard.edu/abs/2001ApJ...563...80S/abstract)
        Eq. 5. and summing over all the bins
    """
    fname = root_dir/'resources/Stern_logRlogN.txt'

    bins, hist_obs, err_obs = read_logRlogN(fname)
    delta_bins = np.log10(bins[1:]/bins[:-1])
    all_sky_glob_rate = np.sum(10**(hist_obs)*delta_bins)

    log.info(''.join([f"Global LGRB rate from Stern constraint: {all_sky_glob_rate:.2f} ",
             f"GRB/yr in 4 pi with peak flux in [50-300 keV] above {bins.min()} ph/s/cm2"]))

    return all_sky_glob_rate


def efficiency_correction_Stern(pflx, c_e0=0.097, nu=2.34, norm=0.7):
    """
        The efficiency function of BATSE for detecting GRBs as a function of peak flux, derived by Stern+01
        c_e0 is in [counts/s/cm2]
        pflx is in [ph/s/cm2]
    """
    c_e = pflx * 0.75  # the conversion factor from counts to pflx comes from the Stern+01 paper as well, figure 7.
    return norm * (1.0 - np.exp(-(c_e/c_e0)**2))**nu


def k_correction(z, photon_index=1.5, zm=2.):
    return ((1.+z)/(1+zm))**(photon_index-2)


def str2bool(s):
    if s.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif s.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise TypeError('Boolean value expected')
