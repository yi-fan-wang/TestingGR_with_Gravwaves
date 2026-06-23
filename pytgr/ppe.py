import re
import numpy

MTSUN_SI = 4.925490947e-6
_PPE_BETA_RE = re.compile(r"^ppebeta(-?\d+)$")

def ppe_beta_exponent(pn_index):
    """Return the ppE phase exponent for a PN-indexed beta parameter."""
    return (int(pn_index) - 5.0) / 3.0

def _extract_ppe_beta_terms(kwds):
    return _extract_indexed_terms(kwds, _PPE_BETA_RE)

def _extract_indexed_terms(kwds, pattern):
    terms = []
    for key in list(kwds):
        match = pattern.match(key)
        if match is None:
            continue
        value = kwds.pop(key)
        if value is None or value == 0:
            continue
        terms.append((int(match.group(1)), value))
    return sorted(terms)

def _ppe_phase_shift(
    frequencies,
    ppe_beta_terms=None,
    chirp_mass_seconds=None,
):
    frequencies = numpy.asarray(frequencies, dtype=float)
    phase = numpy.zeros_like(frequencies, dtype=float)
    has_phase = False

    if ppe_beta_terms:
        if chirp_mass_seconds is None:
            raise ValueError("chirp_mass_seconds is required when ppe_beta_terms are non-zero.")
        u = numpy.pi * chirp_mass_seconds * frequencies
        for pn_index, beta in ppe_beta_terms:
            phase += beta * u**ppe_beta_exponent(pn_index) # ppE phase contribution
        has_phase = True

    if not has_phase:
        return None
    return phase

def gen_ppe_waveform(**kwds):
    '''
    Parameterized Post Einstein

    Parameters
    ----------
    kwds: dict
        Supports original ppE-like phase coefficients ``ppebetaN``. The integer
        suffix ``N`` follows the PN phase index, so ``ppebeta2`` has exponent
        ``b = -1`` and contributes ``ppebeta2 / u`` with
        ``u = pi * detector_chirp_mass * f``.

    Returns
    -------
    hp: pycbc.types.FrequencySeries
        Plus polarization time series
    hc: pycbc.types.FrequencySeries
        Cross polarization time series
    '''
    from pycbc.waveform import get_fd_waveform
    from ._utils import pop_base_gr_approximant, strip_plugin_approximant

    # sanity checks
    base_gr_approximant = pop_base_gr_approximant(kwds)
    ppe_beta_terms = _extract_ppe_beta_terms(kwds)
    chirp_mass_seconds = None
    if ppe_beta_terms:
        from pycbc import conversions

        chirp_mass_seconds = conversions.mchirp_from_mass1_mass2(
            kwds["mass1"], kwds["mass2"]
        )
        chirp_mass_seconds *= MTSUN_SI

    # Generate GR waveforms
    strip_plugin_approximant(kwds)
    hp, hc = get_fd_waveform(approximant=base_gr_approximant, **kwds)

    # Add PPE parameters
    # Slicing with index 1 avoids the zero-frequency bin.
    phase = _ppe_phase_shift(
        hp.sample_frequencies[1:],
        ppe_beta_terms=ppe_beta_terms,
        chirp_mass_seconds=chirp_mass_seconds,
    )
    if phase is not None:
        hp[1:] *= numpy.exp(1j * phase)
        hc[1:] *= numpy.exp(1j * phase)
    return hp, hc
