
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
    from pycbc.types import FrequencySeries
    import lal, lalsimulation as lalsim

    # sanity checks
    if kwds['baseapprox'] is None:
        raise ValueError("A base waveform approximant is required.")
    # Generate GR waveforms
    if 'approximant' in kwds:
        kwds.pop("approximant")
    hp, hc = get_fd_waveform(approximant=kwds['baseapprox'], **kwds)

    nonGRdict = lal.CreateDict()
    if kwds['ftadchi2'] != None:
        lalsim.SimInspiralWaveformParamsInsertNonGRDChi2(nonGRdict,float(kwds['ftadchi2']))
    
    hplal = hp.lal()
    hclal = hc.lal()

    lalsim.SimInspiralTestingGRCorrections(hplal,
                                       2,2,
                                       par['mass1']*lal.MSUN_SI,
                                       par['mass2']*lal.MSUN_SI, 
                                       par['spin1z'], 
                                       par['spin2z'], 
                                       par['f_lower'],#f_start, 
                                       par['f_lower'],#f_ref, 
                                       0.39, #https://arxiv.org/pdf/2109.06988.pdf 
                                       1,
                                       nonGRdict)

    lalsim.SimInspiralTestingGRCorrections(hclal,
                                       2,2,
                                       par['mass1']*lal.MSUN_SI,
                                       par['mass2']*lal.MSUN_SI, 
                                       par['spin1z'], 
                                       par['spin2z'], 
                                       par['f_lower'],#f_start, 
                                       par['f_lower'],#f_ref, 
                                       0.39, #https://arxiv.org/pdf/2109.06988.pdf 
                                       1,
                                       nonGRdict)
    
    # Build the FrequencySeries format
    hpfta = FrequencySeries(hplal.data.data, dtype=hplus.dtype,delta_f=hplal.deltaF, epoch=hplalsim.epoch)
    hcfta = FrequencySeries(hclal.data.data, dtype=hcross.dtype,delta_f=hplal.deltaF, epoch=hclalsim.epoch)
    return hp_mg, hc_mg

def length_in_time(**kwds):
    from pycbc.waveform.waveform import get_waveform_filter_length_in_time
    return get_waveform_filter_length_in_time(**kwds)
