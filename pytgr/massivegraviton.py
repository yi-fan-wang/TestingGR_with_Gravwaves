
def gen_waveform(**kwds):
    '''
    Generate waveform with massive graviton correction, assuming
    no effects are introduced to binary dynamics. Refer to 
    https://journals.aps.org/prd/pdf/10.1103/PhysRevD.57.2061
    "Bounding the mass of the graviton using gravitational-wave 
    observations of inspiralling compact binaries" by Clifford Will
    for more details

    Parameters
    ----------
    kwds: dict
        The parameters defining the waveform to generator. In particular,
    one should provide "lambda_g", the Compton wavelength for massive
    graviton, and 'baseapprox', the based GR waveform, on top of which
    correction from massive graviton will be added.

    Returns
    -------
    hp: pycbc.types.FrequencySeries
        Plus polarization time series
    hc: pycbc.types.FrequencySeries
        Cross polarization time series
    '''
    from pycbc.waveform import get_fd_waveform
    from pycbc import cosmology, pnutils
    import scipy.constants

    # sanity checks
    if kwds['baseapprox'] is None:
        raise ValueError("A base waveform approximant is required.")
    if kwds['lambda_g'] is None:
        raise ValueError("The lambda_g is required for waveform with massive graviton.")

    # Generate GR waveforms
    if 'approximant' in kwds:
        kwds.pop("approximant")
    hp, hc = get_fd_waveform(approximant=kwds['baseapprox'], **kwds)

    # start to comput massive graviton correction terms,
    # phi = - pi * D * c / lambda_g^2 / (1+z) / f
    z = cosmology.redshift(kwds['distance'])
    D = pnutils.megaparsecs_to_meters(float(kwds['distance']))
    c = scipy.constants.speed_of_light
    lambda_g = kwds['lambda_g']
    phi = - numpy.pi * D * c / lambda_g / lambda_g / (1+z) / hp.sample_frequencies

    hp_mg = hp * phi
    hc_mg = hc * phi

    return hp_mg, hc_mg

def length_in_time(**kwds):
    from pycbc.waveform.waveform import get_waveform_filter_length_in_time
    return get_waveform_filter_length_in_time(**kwds)
