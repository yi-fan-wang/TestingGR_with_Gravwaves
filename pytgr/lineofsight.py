import numpy

def gen_waveform(**kwds):
    from pycbc.waveform import get_fd_waveform
    from pycbc import cosmology
    import lal.LAL_MTSUN_SI
    from pycbc import conversions

    if 'approximant' in kwds:
        kwds.pop("approximant")
    if kwds['baseapprox'] is None:
        raise ValueError("A base waveform approximant is required.")
    if kwds['lsc_a'] is None:
        raise ValueError("A line of sight acceleration parameter is required!")

    hp, hc = get_fd_waveform(approximant=kwds['baseapprox'], **kwds)
    fseries = hp.sample_frequencies 

    eta = conversions.eta_from_mass1_mass2(kwds['mass1'],kwds['mass2'])
    m_total = kwds['mass1'] + kwds['mass2']
    v_f = np.pi * lal.LAL_MTSUN_SI * m_total * fseries
    delta_phi = 25 / 65536 / eta / eta * lal.LAL_MTSUN_SI * m_total * kwds['lsc_a'] * vf ** (-13)

    hp_lsa = hp * delta_phi
    hc_lsa = hc * delta_phi

    return hp_lsa, hc_lsa
