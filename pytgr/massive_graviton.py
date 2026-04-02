import numpy as np
import scipy.constants
from scipy.integrate import quad
import pycbc.cosmology

def mg_to_lambda_g(mg):
    '''
    Convert graviton mass to Compton wavelength.

    Parameters
    ----------
    mg : float
        Graviton mass in eV/c^2

    Returns
    -------
    lambda_g : float
        Compton wavelength of graviton in km
    '''
    h_eV_s, _, _ = scipy.constants.physical_constants['Planck constant in eV/Hz']
    return h_eV_s / (mg / scipy.constants.c) / 1000 # convert from m to km

def lambda_g_to_mg(lambda_g):
    '''
    Convert graviton Compton wavelength to mass.
    
    Parameters
    ----------
    lambda_g : float
        Compton wavelength of graviton in km
    
    Returns
    -------
    mg : float
        Graviton mass in eV/c^2
    '''
    h_eV_s, _, _ = scipy.constants.physical_constants['Planck constant in eV/Hz']
    lambda_g_m = lambda_g * 1000 # convert from km to m
    return h_eV_s / (lambda_g_m / scipy.constants.c)

def effective_distance(luminosity_distance):
    '''
    Calculate effective distance D_eff for massive graviton phase correction.
    D_eff = (1+z) * c / H0 * integral_0^z dz' / ((1+z')^2 * E(z'))

    Parameters
    ----------
    luminosity_distance : float
        Luminosity distance to the source in Mpc
    
    Returns
    -------
    deff : float
        Effective distance in Mpc
    '''
    Omega_M = 0.315
    Omega_Lambda = 0.685
    H0 = 67.4  # km/s/Mpc
    H0_m = H0 * 1000  # convert to m/s/Mpc
    z = pycbc.cosmology.redshift(luminosity_distance)

    def integrand(zp):
        E = np.sqrt(Omega_M * (1+zp)**3 + Omega_Lambda)
        return 1.0 / ((1+zp)**2 * E)

    integral, _ = quad(integrand, 0, z)
    deff = scipy.constants.c * (1 + z) / H0_m * integral # in Mpc
    return deff

def mg_phase_correction(mg, distance, frequencies):
    '''
    Calculate the massive graviton phase correction for given frequencies.
    delta_phi = -pi * D_eff * c^3 * m_g^2 / h^2 / (1 + z) / f

    Parameters
    ----------
    mg : float
        Graviton mass in eV/c^2
    distance : float
        Luminosity distance to the source in Mpc
    frequencies : array-like
        Array of frequencies in Hz

    Returns
    -------
    delta_phi : array-like
        Phase correction for each frequency in radians
    '''
    z = pycbc.cosmology.redshift(distance)
    deff = effective_distance(distance) * 1e6 * scipy.constants.parsec # in unit of m
    c = scipy.constants.speed_of_light
    h_eV_s, _, _ = scipy.constants.physical_constants['Planck constant in eV/Hz']
    delta_phi = - np.pi * deff * mg**2 / c / h_eV_s**2 / (1 + z) / frequencies
    return delta_phi

def gen_mg_waveform(**kwds):
    '''
    Generate frequency-domain waveform with massive graviton phase correction.
    Based on Will (1997) PRD 57, 2061. Assumes no modification to binary dynamics.

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
    delta_phi = mg_phase_correction(kwds['mg'], kwds['distance'], hp.sample_frequencies[1:])

    # slicing with index 1 to avoid dividing zero frequency
    hp[1:] *= np.exp(1j * delta_phi)
    hc[1:] *= np.exp(1j * delta_phi)
    
    return hp, hc