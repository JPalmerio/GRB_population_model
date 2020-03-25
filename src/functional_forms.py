import numpy as np


def Schechter_log(logL, logLbreak, slope):
    """
        Returns the unnormalized Schechter function
        Expects Lum arguments to be in log scale
    """
    x = 10.**(logL - logLbreak)
    Sch = x**(1.-slope) * np.exp(-x)
    return Sch


def BPL_lum(logL, logLbreak, slopeL, slopeH):
    """
        Returns the unnormalized broken power law function
        Expects Lum arguments to be in log scale
    """
    x = 10.**(logL - logLbreak)
    BPL_func = np.where(x <= 1, x**(1.-slopeL), x**(1.-slopeH))
    return BPL_func


def SH(z, a=2.37, b=1.8, zm=2, nu=0.178, IMF_norm=0.007422):
    """
        Springel-Hernquist+03 functional form for the cosmic SFR.
        Default are given the values of Vangioni+15.
        Returns an event rate in units of yr-1 Mpc-3
        Note : nu is in units of Msun/yr/Mpc3 and IMF_norm in units of Msun-1
    """
    return IMF_norm * nu * a * np.exp(b*(z-zm)) / ((a-b) + b*np.exp(a*(z-zm)))


def BExp(z, a=1.1, b=-0.57, zm=1.9, SFR_norm=0.02744, IMF_norm=0.007422):
    """
        GRB rate as parametrized by a broken exponential function.
        Default values are chosen as best fit from SFR of Vangioni+15
        If you leave the default SFR_norm and IMF_norm, the result
        will be in units of yr-1 Mpc-3.
        IMF_norm is in units of M_sun-1 and converts the CSFRD (in
        units of M_sun yr-1 Mpc-3) to a core-collapse rate density (in
        units of yr-1 Mpc-3).
        SFR_norm is adjusted on the functional form of Springel-Hernquist
        (SH) with the parameter values of Vangioni+15.
    """
    rate = np.where(z <= zm, np.exp(a*z), np.exp(b*z) * np.exp((a-b)*zm))
    return SFR_norm*IMF_norm*rate


def BPL_z(z, a=2.07, b=-1.36, zm=3.11, norm=1.3e-9):
    """
        Returns the unnormalized broken power law function
        Expects Lum arguments to be in log scale
        The default values are from Wanderman & Piran 2010
        The norm is given by them as well, converted to units of
        yr-1 Mpc-3 (carefull this is the LGRB rate not the core-
        collapse rate, you need to multiply by the efficiency and
        the average opening angle.)
    """
    rate = np.where(z <= zm, (1.+z)**a, (1.+z)**b * (1.+zm)**(a-b))
    return norm*rate


def MD(z, gamma_0=0.0204, gamma_1=1.8069, gamma_2=3.1724, gamma_3=7.2690):
    """
        Fit of the GRB rate from Pescalli et al. 2016 based on the
        Star Formation Rate Density from Madau & Dickinson 2014.
        Normalized to its maximum (arbitrary units)
    """
    return gamma_0 * (1.+z)**gamma_1 / (1. + ((1.+z)/gamma_2)**gamma_3)


def D06(z, a, b, c, d, IMF_norm=0.007422):
    """
        SFR of Daigne+06
        Default values for [SFR1, SFR2, SFR3] are:
        a = [0.320, 0.196, 0.175]
        b = [3.30, 4.0, 3.67]
        c = [3.52, 4.0, 3.58]
        d = [23.6, 14.6, 12.6]
    """
    return IMF_norm * a*np.exp(b*z) / (d + np.exp(c*z))
