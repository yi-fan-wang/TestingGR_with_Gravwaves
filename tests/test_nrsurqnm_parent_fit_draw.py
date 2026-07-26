import unittest

import numpy as np

try:
    from pycbc.waveform import get_td_waveform
    from tgr.nrsurqnm import FIT_TSTART_MIN
except ImportError:
    get_td_waveform = None
    FIT_TSTART_MIN = None

COMMON = dict(
    approximant='NRSur7dq4_remove_qqnm',
    mass1=35.0, mass2=30.0,
    spin1x=0.05, spin1y=0.02, spin1z=0.1,
    spin2x=0.01, spin2y=0.03, spin2z=0.05,
    distance=1000.0, inclination=0.5, coa_phase=1.0,
    delta_t=1.0 / 4096, f_lower=20.0, f_ref=20.0,
    mode22='220 221 222 223',
    mode22_omitted='224',
    mode_quadratic='220220 220221 221221 220222 221222 222222 '
                   '220223 221223 222223 223223',
)
# below FIT_TSTART_MIN so the grid-fit branch (where parent_fit_draw acts) runs
TOFFSET_EARLY = 0.001


def _nrsur_available():
    if get_td_waveform is None:
        return False
    try:
        get_td_waveform(toffset=FIT_TSTART_MIN, **COMMON)
    except Exception:
        return False
    return True


@unittest.skipUnless(_nrsur_available(),
                     "pycbc/lal or NRSur7dq4 data not available")
class ParentFitDrawTests(unittest.TestCase):
    def test_mean_is_deterministic(self):
        hp1, _ = get_td_waveform(toffset=TOFFSET_EARLY,
                                 parent_fit_draw='mean', **COMMON)
        hp2, _ = get_td_waveform(toffset=TOFFSET_EARLY,
                                 parent_fit_draw='mean', **COMMON)
        self.assertTrue(np.array_equal(hp1.numpy(), hp2.numpy()))

    def test_gaussian_default_is_stochastic(self):
        # default draw ('gaussian', no seed) changes between calls
        hp1, _ = get_td_waveform(toffset=TOFFSET_EARLY, **COMMON)
        hp2, _ = get_td_waveform(toffset=TOFFSET_EARLY, **COMMON)
        self.assertFalse(np.array_equal(hp1.numpy(), hp2.numpy()))

    def test_gaussian_seeded_is_reproducible(self):
        hp1, _ = get_td_waveform(toffset=TOFFSET_EARLY, seed=42, **COMMON)
        hp2, _ = get_td_waveform(toffset=TOFFSET_EARLY, seed=42, **COMMON)
        self.assertTrue(np.array_equal(hp1.numpy(), hp2.numpy()))

    def test_mean_differs_from_gaussian_draw(self):
        hp_mean, _ = get_td_waveform(toffset=TOFFSET_EARLY,
                                     parent_fit_draw='mean', **COMMON)
        hp_draw, _ = get_td_waveform(toffset=TOFFSET_EARLY, seed=42, **COMMON)
        self.assertFalse(np.array_equal(hp_mean.numpy(), hp_draw.numpy()))

    def test_ignored_above_fit_tstart_min(self):
        # single direct fit at t0: parent_fit_draw must have no effect
        hp1, _ = get_td_waveform(toffset=FIT_TSTART_MIN,
                                 parent_fit_draw='mean', **COMMON)
        hp2, _ = get_td_waveform(toffset=FIT_TSTART_MIN,
                                 parent_fit_draw='gaussian', **COMMON)
        self.assertTrue(np.array_equal(hp1.numpy(), hp2.numpy()))

    def test_invalid_value_raises(self):
        with self.assertRaises(ValueError):
            get_td_waveform(toffset=TOFFSET_EARLY,
                            parent_fit_draw='median', **COMMON)


if __name__ == '__main__':
    unittest.main()
