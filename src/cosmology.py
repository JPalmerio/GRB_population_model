import numpy as np
from constants import cLight


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


def create_cosmology(OmegaM=0.27, OmegaL=0.73, h=0.71, zmax=20, z_step=0.001, verbose=True):
    """
        Create your own cosmology with OmegaM, OmegaL and H0 = 100*h.
        step is the interval of the redshift array.
    """
    if verbose:
        print("Creating cosmology with: OmegaM = {}, OmegaL = {}, H0 = {}".format(OmegaM, OmegaL, 100.*h))
    cosmo = {}
    z = np.arange(0, zmax, z_step)
    N = len(z)
    cdivH0 = cLight/(10.**7*h)
    # Luminosity distance
    D_L = np.zeros(N)
    for i in range(1,N):
        D_L[i] = (1+z[i])/(1+z[i-1])*D_L[i-1] + (1+z[i])*cdivH0 * 0.5*z_step\
            * (1./E(z[i], OmegaM=OmegaM, OmegaL=OmegaL) + 1./E(z[i-1], OmegaM=OmegaM, OmegaL=OmegaL))

    # Differential comoving volume element
    dVdz = np.zeros(N)
    for i in range(1,N):
        dVdz[i] = 4.*np.pi * cdivH0 * D_L[i]**2/(1.+z[i])**2 * 1./E(z[i], OmegaM=OmegaM, OmegaL=OmegaL)

    # Volume of the Universe at redshift z
    Vz = np.zeros(N)
    for i in range(1,N):
        Vz[i] = Vz[i-1] + 0.5*z_step * (dVdz[i]/(1.+z[i]) + dVdz[i-1]/(1.+z[i-1]))

    cosmo['redshift'] = z
    cosmo['D_L'] = D_L
    cosmo['dVdz'] = dVdz
    cosmo['Vz'] = Vz
    return cosmo


def dVdz(z, cosmo):
    """ Returns the differential comoving volume element at a given z """
    z_ind = cosmo['redshift'].searchsorted(z)
    dVdz = cosmo['dVdz'][z_ind]
    return dVdz


def Lum_dist(z, cosmo):
    """ Returns Luminosity Distance found in D_L.dat file (in Mpc) """

    z_ind = cosmo['redshift'].searchsorted(z)
    Lum_dist = cosmo['D_L'][z_ind]

    return Lum_dist


def E(z, OmegaM, OmegaL):
    """ True for OmegaK = 0 """
    E = np.sqrt(OmegaM*(1+z)**3 + OmegaL)
    return E
