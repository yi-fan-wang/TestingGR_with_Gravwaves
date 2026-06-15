
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
    from ._utils import pop_base_gr_approximant, strip_plugin_approximant
    import numpy

    # sanity checks
    base_gr_approximant = pop_base_gr_approximant(kwds)
    # Generate GR waveforms
    strip_plugin_approximant(kwds)
    hp, hc = get_fd_waveform(approximant=base_gr_approximant, **kwds)

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
