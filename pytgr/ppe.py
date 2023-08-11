
def gen_waveform(**kwds):
    '''
    Parameterized Post Einstein

    Parameters
    ----------
    kwds: dict
        Only support ppedchi2 atm.

    Returns
    -------
    hp: pycbc.types.FrequencySeries
        Plus polarization time series
    hc: pycbc.types.FrequencySeries
        Cross polarization time series
    '''
    from pycbc.waveform import get_fd_waveform
    import numpy

    # sanity checks
    if kwds['baseapprox'] is None:
        raise ValueError("A base waveform approximant is required.")
    # Generate GR waveforms
    if 'approximant' in kwds:
        kwds.pop("approximant")
    hp, hc = get_fd_waveform(approximant=kwds['baseapprox'], **kwds)

    # Add PPE parameters
    if kwds['ppedchi2'] != None:
        if kwds['ppedchi2'] == 0:
            return hp, hc
        else:
            phi = kwds['ppedchi2'] / hp.sample_frequencies[1:]
            # slicing with index 1 to avoid dividing zero frequency

    hp[1:] *= numpy.exp(-1j*phi)
    hc[1:] *= numpy.exp(-1j*phi)
    return hp, hc
