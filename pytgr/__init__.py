def length_in_time(**kwds):
    from pycbc.waveform.waveform import get_waveform_filter_length_in_time
    if "approximant" in kwds:
        kwds.pop("approximant")
    return get_waveform_filter_length_in_time(approximant=kwds["baseapprox"], **kwds)
