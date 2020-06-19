import numpy as np
import logging
import scipy.integrate as integrate
from io_grb_pop import read_column

log = logging.getLogger(__name__)

try:
    import f90_functions as f90f
    f90 = True
except ImportError:
    f90 = False
    log.error("Could not import f90_functions, f90 set to False")


def init_ECLAIRs(ECLAIRs_dir, Emin, Emax, n_sigma=None, f90=True):
    """
        Function to estimate detection threshold of ECLAIRs using
        realistic noise background.
        Emin and Emax are assumed to be in keV.
        If no sigma was set, the default 6.5 is assumed.
    """

    if n_sigma is None:
        n_sigma = 6.5

    log.debug("==============================================")
    log.debug("===           ECLAIRs instrument           ===")
    log.debug("==============================================")
    log.debug(f"ECLAIRs energy channel: {Emin} to {Emax} keV")
    log.debug(f"ECLAIRs detection level: {n_sigma} sigmas")

    eff_area_E = read_column(ECLAIRs_dir/'ECLAIRs_rf_eff.txt', 0)
    eff_area_A = read_column(ECLAIRs_dir/'ECLAIRs_rf_eff.txt', 1)
    bkg_E = read_column(ECLAIRs_dir/'ECLAIRs_rf_bkg.txt', 0)
    bkg_B = read_column(ECLAIRs_dir/'ECLAIRs_rf_bkg.txt', 1)
    offax_corr = np.genfromtxt(ECLAIRs_dir/'ECLAIRs_rf_offaxis.txt')
    omega_ECLAIRs = np.genfromtxt(ECLAIRs_dir/'ECLAIRs_rf_omega.txt')

    Erg_mask = np.where((eff_area_E >= Emin) & (eff_area_E <= Emax))[0]
    if f90:
        bkg_total = f90f.f90f.integrate_1d(x=bkg_E[Erg_mask], y=bkg_B[Erg_mask])
    else:
        bkg_total = integrate.trapz(bkg_B[Erg_mask], bkg_E[Erg_mask], axis=0)

    log.debug(f"Effective area: from {eff_area_A[0]:e} cm2 at {eff_area_E[0]:.2f} keV")
    log.debug(f"Effective area: to   {eff_area_A[-1]:e} cm2 at {eff_area_E[-1]:.2f} keV")
    log.debug(f"Background: from {bkg_B[0]:e} cts/s/keV at {bkg_E[0]:.2f} keV")
    log.debug(f"Background: to   {bkg_B[-1]:e} cts/s/keV at {bkg_E[-1]:.2f} keV")
    log.debug(f"ECLAIRs background: {bkg_total:e} cts/s from {Emin:.2f} keV to {Emax:.2f} keV")

    ECLAIRs_prop = {'eff_area_E': eff_area_E,
                    'eff_area_A': eff_area_A,
                    'bkg_E': bkg_E,
                    'bkg_B': bkg_B,
                    'bkg_total': bkg_total,
                    'offax_corr': offax_corr,
                    'omega_ECLAIRs': omega_ECLAIRs,
                    'omega_ECLAIRs_tot': omega_ECLAIRs.sum(),
                    'n_sigma': n_sigma,
                    'Emin': Emin,
                    'Emax': Emax}

    log.debug("==============================================")

    return ECLAIRs_prop
