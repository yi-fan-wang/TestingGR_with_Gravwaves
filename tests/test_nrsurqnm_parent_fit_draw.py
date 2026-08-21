import unittest

import numpy as np

try:
    import lal
    from pycbc.waveform import get_td_waveform
    from tgr.nrsurqnm import (FIT_TSTART_GRID, FIT_TSTART_MIN, QNMTable,
                              _mode_label_list, _parent_fit_start_grid)
except ImportError:
    lal = None
    get_td_waveform = None
    FIT_TSTART_GRID = None
    FIT_TSTART_MIN = None
    QNMTable = None
    _mode_label_list = None
    _parent_fit_start_grid = None

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

    def test_mf_grid_seeded_is_reproducible_and_changes_waveform(self):
        kw = dict(COMMON, parent_fit_grid_mf='6 7 8 9 10',
                  wls_epsilon_floor=1e-22)
        hp1, _ = get_td_waveform(toffset=TOFFSET_EARLY, seed=42, **kw)
        hp2, _ = get_td_waveform(toffset=TOFFSET_EARLY, seed=42, **kw)
        hp_old, _ = get_td_waveform(toffset=TOFFSET_EARLY, seed=42, **COMMON)
        self.assertTrue(np.array_equal(hp1.numpy(), hp2.numpy()))
        self.assertFalse(np.array_equal(hp1.numpy(), hp_old.numpy()))

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

    def test_direct_is_deterministic(self):
        hp1, _ = get_td_waveform(toffset=TOFFSET_EARLY,
                                 parent_fit_draw='direct', **COMMON)
        hp2, _ = get_td_waveform(toffset=TOFFSET_EARLY,
                                 parent_fit_draw='direct', **COMMON)
        self.assertTrue(np.array_equal(hp1.numpy(), hp2.numpy()))

    def test_direct_differs_from_mean_below_fit_tstart_min(self):
        hp_direct, _ = get_td_waveform(toffset=TOFFSET_EARLY,
                                       parent_fit_draw='direct', **COMMON)
        hp_mean, _ = get_td_waveform(toffset=TOFFSET_EARLY,
                                     parent_fit_draw='mean', **COMMON)
        self.assertFalse(np.array_equal(hp_direct.numpy(), hp_mean.numpy()))

    def test_direct_equals_single_fit_path(self):
        # 'direct' below FIT_TSTART_MIN must reproduce the ordinary
        # single-fit-at-t0 branch, i.e. what one gets by lowering
        # FIT_TSTART_MIN under the same start time
        import tgr.nrsurqnm as nr
        hp_direct, _ = get_td_waveform(toffset=TOFFSET_EARLY,
                                       parent_fit_draw='direct', **COMMON)
        saved = nr.FIT_TSTART_MIN
        try:
            nr.FIT_TSTART_MIN = 0.5 * TOFFSET_EARLY
            hp_low, _ = get_td_waveform(toffset=TOFFSET_EARLY,
                                        parent_fit_draw='mean', **COMMON)
        finally:
            nr.FIT_TSTART_MIN = saved
        self.assertTrue(np.array_equal(hp_direct.numpy(), hp_low.numpy()))

    def test_direct_ignored_above_fit_tstart_min(self):
        hp1, _ = get_td_waveform(toffset=FIT_TSTART_MIN,
                                 parent_fit_draw='direct', **COMMON)
        hp2, _ = get_td_waveform(toffset=FIT_TSTART_MIN,
                                 parent_fit_draw='mean', **COMMON)
        self.assertTrue(np.array_equal(hp1.numpy(), hp2.numpy()))


@unittest.skipUnless(_mode_label_list is not None, "tgr not available")
class ModeLabelCoercionTests(unittest.TestCase):
    # inference config files convert numeric-looking static params to
    # float, so a single mode label arrives as e.g. 224.0
    def test_label_list_forms(self):
        self.assertEqual(_mode_label_list('220 221 222'), ['220', '221', '222'])
        self.assertEqual(_mode_label_list('224'), ['224'])
        self.assertEqual(_mode_label_list(224.0), ['224'])
        self.assertEqual(_mode_label_list(224), ['224'])
        self.assertEqual(_mode_label_list(None), [])
        with self.assertRaises(ValueError):
            _mode_label_list(224.5)

    @unittest.skipUnless(_nrsur_available(),
                         "pycbc/lal or NRSur7dq4 data not available")
    def test_float_omitted_matches_string(self):
        kw = dict(COMMON)
        kw['mode22_omitted'] = 224.0
        hp_float, _ = get_td_waveform(toffset=FIT_TSTART_MIN, **kw)
        hp_str, _ = get_td_waveform(toffset=FIT_TSTART_MIN, **COMMON)
        self.assertTrue(np.array_equal(hp_float.numpy(), hp_str.numpy()))

    @unittest.skipUnless(_nrsur_available(),
                         "pycbc/lal or NRSur7dq4 data not available")
    def test_float_seed_matches_int_seed(self):
        hp_f, _ = get_td_waveform(toffset=TOFFSET_EARLY, seed=42.0, **COMMON)
        hp_i, _ = get_td_waveform(toffset=TOFFSET_EARLY, seed=42, **COMMON)
        self.assertTrue(np.array_equal(hp_f.numpy(), hp_i.numpy()))


@unittest.skipUnless(_parent_fit_start_grid is not None, "tgr not available")
class ParentFitStartGridTests(unittest.TestCase):
    def setUp(self):
        self.qnm_par = QNMTable(final_mass=100.0, final_spin=0.7,
                                freq={}, tau={})

    def test_default_retains_historical_seconds_grid(self):
        np.testing.assert_array_equal(
            _parent_fit_start_grid(self.qnm_par), FIT_TSTART_GRID)

    def test_mf_grid_uses_sample_final_mass(self):
        expected = np.arange(6.0, 11.0) * 100.0 * lal.MTSUN_SI
        np.testing.assert_allclose(
            _parent_fit_start_grid(self.qnm_par, '6 7 8 9 10'), expected,
            rtol=0, atol=0)

    def test_mf_grid_scales_with_final_mass(self):
        other = QNMTable(final_mass=80.0, final_spin=0.7, freq={}, tau={})
        grid_100 = _parent_fit_start_grid(self.qnm_par, [6, 7, 8, 9, 10])
        grid_80 = _parent_fit_start_grid(other, [6, 7, 8, 9, 10])
        np.testing.assert_allclose(grid_80 / grid_100, 0.8,
                                   rtol=1e-15, atol=0)

    def test_invalid_mf_grids_raise(self):
        for grid in ('6', '6 6 8', '6 5', '0 6', '6 nan'):
            with self.subTest(grid=grid), self.assertRaises(ValueError):
                _parent_fit_start_grid(self.qnm_par, grid)


if __name__ == '__main__':
    unittest.main()
