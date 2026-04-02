import numpy as np
import scipy.constants
from scipy.integrate import quad
import pycbc.cosmology

def mg_to_lambda_g(mg):
    return scipy.constants.h / (mg * scipy.constants.c)

def lambda_g_to_mg(lambda_g):
    return scipy.constants.h / (lambda_g * scipy.constants.c)

def effective_distance(luminosity_distance):
    Omega_M = 0.315
    Omega_Lambda = 0.685
    H0 = 67.4  # km/s/Mpc
    H0_m = H0 * 1000  # convert to m/s/Mpc
    z = pycbc.cosmology.redshift(luminosity_distance)

    def integrand(zp):
        E = np.sqrt(Omega_M * (1+zp)**3 + Omega_Lambda)
        return 1.0 / ((1+zp)**2 * E)

    integral, _ = quad(integrand, 0, z)
    D = scipy.constants.c * (1 + z) / H0_m * integral # in Mpc
    return D

def gen_mg_waveform(**kwds):
    '''
    Generate frequency-domain waveform with massive graviton phase correction.

    Based on Will (1997) PRD 57, 2061. Assumes no modification to binary dynamics.
    Phase correction: delta_phi = -pi * D_c * c / lambda_g^2 / f

    Parameters
    ----------
    kwds : dict
        Must contain:
            baseapprox : str        GR approximant name
            mg : float              graviton mass in eV/c^2
            distance : float        luminosity distance in Mpc
        Other parameters (masses, spins, etc.) passed to get_fd_waveform.

    Returns
    -------
    hp: pycbc.types.FrequencySeries
        Plus polarization time series
    hc: pycbc.types.FrequencySeries
        Cross polarization time series
    '''
    from pycbc.waveform import get_fd_waveform

    # sanity checks
    for key in ['base_gr_approximant', 'mg', 'distance']:
        if key not in kwds or kwds[key] is None:
            raise ValueError(f"Missing required argument: {key}")

    # Generate GR waveforms
    if 'approximant' in kwds:
        kwds.pop("approximant")
    hp, hc = get_fd_waveform(approximant=kwds['base_gr_approximant'], **kwds)

    # Apply massive graviton phase correction
    z = pycbc.cosmology.redshift(kwds['distance'])
    deff = effective_distance(kwds['distance'])
    mg = kwds['mg'] # in unit of eV/c^2
    c = scipy.constants.speed_of_light
    h_eV_s, _, _ = scipy.constants.physical_constants['Planck constant in eV/Hz']
    delta_phi = - np.pi * deff * mg**2 * c / h_eV_s**2 / (1 + z) / hp.sample_frequencies[1:]

    # slicing with index 1 to avoid dividing zero frequency
    hp[1:] *= np.exp(1j * delta_phi)
    hc[1:] *= np.exp(1j * delta_phi)

    return hp, hc