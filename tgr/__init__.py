def length_in_time(**kwds):
    from pycbc.waveform.waveform import get_waveform_filter_length_in_time
    from ._utils import pop_base_gr_approximant, strip_plugin_approximant

    strip_plugin_approximant(kwds)
    base_gr_approximant = pop_base_gr_approximant(kwds)
    return get_waveform_filter_length_in_time(approximant=base_gr_approximant, **kwds)
