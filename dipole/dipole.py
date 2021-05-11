from astropy.constants import G,c,M_sun
s_g = (G/c**3)

def IMRPhenomdipole(B=0.0, **kwds):
    from pycbc.waveform import get_fd_waveform
    import numpy as np
    import pycbc.conversions

    if 'approximant' in kwds:
        kwds.pop("approximant")
    hp, hc = get_fd_waveform(approximant="IMRPhenomXPHM", **kwds)

    eta = pycbc.conversions.eta_from_mass1_mass2(kwds['mass1'],kwds['mass2'])
    M_chirp = pycbc.conversions.mchirp_from_mass1_mass2(kwds['mass1'],kwds['mass2'])*M_sun.value
    M = pycbc.conversions.mtotal_from_mass1_mass2(kwds['mass1'],kwds['mass2'])*M_sun.value

    fISCO = 1/(6**(3/2)*pi*M*s_g)

    fsampling = hp.sample_frequencies

    beta = -3/224*eta**(2/5)*B
    b = -7/3
    dipole = np.exp(1j*beta*(pi*M_chirp*fsampling*s_g)**b)
    dipole[(fsampling>fISCO)] = np.ones_like(fsampling[(fsampling>fISCO)])

    hpd, hcd = hp.data*dipole, hc.data*dipole

    #derive_B = (-1j) * (3/(224 * eta)) * (pi*M*fsampling*s_g)**(-7/3) * hcd

    return hpd,hcd