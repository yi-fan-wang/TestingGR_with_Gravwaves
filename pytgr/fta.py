
def gen_waveform(**kwds):
    '''
    Generate waveform with FTA (flexible theory agnostic) correction.
    Described in https://arxiv.org/pdf/1811.00364.pdf "Tests of 
    General Relativity with GW170817"

    Parameters
    ----------
    kwds: dict
        Only support dchi2 atm.

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
    # add FTA correction
    lalsim.SimInspiralTestingGRCorrections(hplal,
                                       2,2, #only support (2,2) mode
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
    
    # build FrequencySeries format
    hpfta = FrequencySeries(hplal.data.data, dtype=hplus.data.data.dtype,
                            delta_f=hplal.deltaF, epoch=hplal.epoch)
    hcfta = FrequencySeries(hclal.data.data, dtype=hcross.data.data.dtype,
                            delta_f=hplal.deltaF, epoch=hclal.epoch)
    return hpfta, hcfta
