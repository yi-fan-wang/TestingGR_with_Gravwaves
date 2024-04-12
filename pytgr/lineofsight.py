def gen_waveform(**kwds):
    from pycbc.waveform import get_fd_waveform
    import lal, numpy
    from pycbc import conversions

    if 'approximant' in kwds:
        kwds.pop("approximant")
    if kwds['baseapprox'] is None:
        raise ValueError("A base waveform approximant is required.")
    if kwds['lsa_a'] is None:
        raise ValueError("A line of sight acceleration parameter is required!")

    hp, hc = get_fd_waveform(approximant=kwds['baseapprox'], **kwds)
    fseries = hp.sample_frequencies 

    eta = conversions.eta_from_mass1_mass2(kwds['mass1'],kwds['mass2'])
    m_total = kwds['mass1'] + kwds['mass2']
    v_f = lal.PI * lal.MTSUN_SI * m_total * fseries[1:]
    delta_phi = 25 / 65536 / eta / eta * lal.MTSUN_SI * m_total * kwds['lsa_a'] / lal.C_SI * v_f ** (-13)

    hp_lsa = hp[1:] * numpy.exp(1j * delta_phi)
    hc_lsa = hc[1:] * numpy.exp(1j * delta_phi)

    return hp_lsa, hc_lsa
