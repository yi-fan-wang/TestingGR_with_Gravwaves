def genwav(dipole_b=0.0, **kwds):
    from pycbc.waveform import get_fd_waveform
    import numpy as np
    import pycbc.conversions as conversions
    import lal
    if 'approximant' in kwds:
        kwds.pop("approximant")
    hp, hc = get_fd_waveform(approximant="IMRPhenomXPHM", **kwds)

    eta = conversions.eta_from_mass1_mass2(kwds['mass1'],kwds['mass2'])
    M_chirp = conversions.mchirp_from_mass1_mass2(kwds['mass1'],kwds['mass2'])
    M = conversions.mtotal_from_mass1_mass2(kwds['mass1'],kwds['mass2'])

    
    kmin = int(kwds['f_lower']/kwds['delta_f'])
    fsampling = hp.sample_frequencies[kmin:]

    beta = -3/224*eta**(2/5)*dipole_b
    dipole = np.exp(1j*beta*(np.pi*M_chirp*fsampling*lal.MTSUN_SI)**(-7/3))
    #dipole[(fsampling>fISCO)] = np.ones_like(fsampling[(fsampling>fISCO)])
    
    hp[kmin:], hc[kmin:] = hp[kmin:]*dipole, hc[kmin:]*dipole

    #derive_B = (-1j) * (3/(224 * eta)) * (pi*M*fsampling*s_g)**(-7/3) * hcd

    return hp,hc


def genwav_seobnrv4_rom(dipole_b=0.0, **kwds):
    from pycbc.waveform import get_fd_waveform
    import numpy as np
    import pycbc.conversions as conversions
    import lal
    if 'approximant' in kwds:
        kwds.pop("approximant")
    hp, hc = get_fd_waveform(approximant="SEOBNRv4_ROM", **kwds)

    eta = conversions.eta_from_mass1_mass2(kwds['mass1'],kwds['mass2'])
    M_chirp = conversions.mchirp_from_mass1_mass2(kwds['mass1'],kwds['mass2'])
    M = conversions.mtotal_from_mass1_mass2(kwds['mass1'],kwds['mass2'])

    
    kmin = int(kwds['f_lower']/kwds['delta_f'])
    fsampling = hp.sample_frequencies[kmin:]

    beta = -3/224*eta**(2/5)*dipole_b
    dipole = np.exp(1j*beta*(np.pi*M_chirp*fsampling*lal.MTSUN_SI)**(-7/3))
    #dipole[(fsampling>fISCO)] = np.ones_like(fsampling[(fsampling>fISCO)])
    
    hp[kmin:], hc[kmin:] = hp[kmin:]*dipole, hc[kmin:]*dipole

    #derive_B = (-1j) * (3/(224 * eta)) * (pi*M*fsampling*s_g)**(-7/3) * hcd

    return hp,hc