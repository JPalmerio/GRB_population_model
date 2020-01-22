import numpy as np


def init_cosmology(cosmo_dir):
    """
        Reads the files created by cosmology.f90 and loads them into a dict
    """

    cosmo = {}
    quant_names = ['redshift', 'D_L', 'dVdz', 'Vz']
    cosmo_quantities = ['precisez', 'D_L', 'dVdz', 'Vz']
    for quant, name in zip(cosmo_quantities, quant_names):
        with open(cosmo_dir/f'Tab{quant}.dat', 'rb') as f:
            cosmo[name] = np.fromfile(f, dtype=np.dtype("i4,100001f8,i4"), count=-1)[0][1]

    return cosmo


def dVdz(z, cosmo):
    """ Returns the differential comoving volume element at a given z"""
    z_ind = cosmo['redshift'].searchsorted(z)
    dVdz = cosmo['dVdz'][z_ind]
    return dVdz


def Lum_dist(z, cosmo):
    """ Returns Luminosity Distance found in D_L.dat file (in Mpc)"""

    z_ind = cosmo['redshift'].searchsorted(z)
    Lum_dist = cosmo['D_L'][z_ind]

    return Lum_dist
