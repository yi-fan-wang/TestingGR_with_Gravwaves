import numpy

def integrand(redshift):
    """
    The integrand:
    (1.0 + z)^parity_beta / sqrt(Omega_m (1+z)^3 + Omega_Lambda)
    """
    omega_m = 0.3075 #pycbc.cosmology.get_cosmology().Om0 # matter density
    omega_l = 0.6910098821161554 #pycbc.cosmology.get_cosmology().Ode0 # dark energy density

    return (1.0+redshift)/ numpy.sqrt(omega_m*(1.0+redshift)**3.0 + omega_l)


def gen_waveform(**kwds):
    from pycbc.waveform import get_fd_waveform
    from pycbc import cosmology
    import lal
    from scipy import integrate

    if 'approximant' in kwds:
        kwds.pop("approximant")
    if kwds['baseapprox'] is None:
        raise ValueError("A base waveform approximant is required.")

    hp, hc = get_fd_waveform(approximant=kwds['baseapprox'], **kwds)
    zz = cosmology.redshift(kwds['distance'])
    intz = integrate.quad(integrand, 0, zz)[0]

    temp =  kwds['parity_mpvinverse'] * intz / 1e9 / lal.QE_SI * (lal.H_SI / 2 / lal.PI) * lal.PI * lal.PI / lal.H0_SI
    expminus = numpy.exp(-1j * temp * hp.sample_frequencies**2)
    expplus = 1 / expminus

    hp_parity = (hp + 1j*hc) * expminus / 2 + (hp - 1j*hc) * expplus / 2
    hc_parity = (hp + 1j*hc) * expminus / 2j - (hp - 1j*hc) * expplus / 2j

    return hp_parity, hc_parity
