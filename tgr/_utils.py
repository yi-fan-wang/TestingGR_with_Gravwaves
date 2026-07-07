def strip_plugin_approximant(kwds):
    """Remove the plugin approximant before calling the GR waveform backend."""
    kwds.pop("approximant", None)


def pop_base_gr_approximant(kwds):
    """Return the GR base approximant configured for plugin waveforms."""
    if "base_gr_approximant" not in kwds or kwds["base_gr_approximant"] is None:
        raise ValueError("A base waveform approximant is required via base_gr_approximant.")
    return kwds.pop("base_gr_approximant")
