"""Tests for remove_fourfour (NOFOURFOUR model) and the free 441
(add441_amp / add441_phi) in NRSur7dq4_remove_qqnm."""
import unittest

import numpy as np

try:
    import lal
    from pycbc.waveform import get_td_waveform, get_td_waveform_modes
    from tgr.nrsurqnm import get_qnm_freqtau
except ImportError:
    lal = None
    get_td_waveform = None
    get_td_waveform_modes = None
    get_qnm_freqtau = None

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
    parent_fit_draw='direct',   # deterministic below FIT_TSTART_MIN
)
TOFFSET = 0.001
AMP, PHI = 0.1, 0.7   # geometric r|A|/M_f and radians


def _nrsur_available():
    if get_td_waveform is None:
        return False
    try:
        get_td_waveform(toffset=TOFFSET, **COMMON)
    except Exception:
        return False
    return True


def _wf(**extra):
    hp, hc = get_td_waveform(toffset=TOFFSET, **dict(COMMON, **extra))
    return hp, hc


def _expected_pol44(series44):
    """(4,4)-sector plus/cross from a complex series, as in the generator."""
    y44 = lal.SpinWeightedSphericalHarmonic(
        COMMON['inclination'], np.pi / 2 - COMMON['coa_phase'], -2, 4, 4)
    y4m4 = lal.SpinWeightedSphericalHarmonic(
        COMMON['inclination'], np.pi / 2 - COMMON['coa_phase'], -2, 4, -4)
    pol = series44 * y44 + np.conj(series44) * y4m4
    return pol.real, -pol.imag


@unittest.skipUnless(_nrsur_available(),
                     "pycbc/lal or NRSur7dq4 data not available")
class Free441Tests(unittest.TestCase):
    def test_defaults_do_not_change_waveform(self):
        # explicit remove_fourfour=0 / add441_amp=0 == parameters absent
        hp0, _ = _wf()
        hp1, _ = _wf(remove_fourfour=0, add441_amp=0.0, add441_phi=PHI)
        self.assertTrue(np.array_equal(hp0.numpy(), hp1.numpy()))

    def test_remove_fourfour_drops_44_sector(self):
        # base = subtract nothing (quadratic_tgr=0): 44 sector kept from t0;
        # no44 = remove_fourfour: difference must be exactly the projected
        # raw-NRSur (4,4) slice
        hp_base, _ = _wf(quadratic_tgr=0.0)
        hp_no44, _ = _wf(remove_fourfour=1)
        # same mode_array as the generator so the time grid is identical
        hlm = get_td_waveform_modes(
            approximant='NRSur7dq4',
            **{k: COMMON[k] for k in ('mass1', 'mass2', 'spin1x', 'spin1y',
                                      'spin1z', 'spin2x', 'spin2y', 'spin2z',
                                      'distance', 'delta_t', 'f_lower',
                                      'f_ref')},
            mode_array=['22', '21', '20', '33', '32', '31', '30',
                        '44', '43', '42', '41', '40'])
        h44 = hlm[(4, 4)][0] + 1j * hlm[(4, 4)][1]
        i0 = int(np.floor(float(TOFFSET - h44.start_time) * h44.sample_rate))
        exp_p = np.zeros(len(h44))
        # the generator's time_slice is end-exclusive: last sample not touched
        exp_p[i0:-1], _ = _expected_pol44(h44.numpy()[i0:-1])
        diff = hp_base.numpy() - hp_no44.numpy()
        scale = np.abs(hp_base.numpy()).max()
        self.assertTrue(np.allclose(diff, exp_p, atol=1e-10 * scale))

    def test_add441_matches_analytic_sinusoid(self):
        hp0, hc0 = _wf(remove_fourfour=1)
        hp1, hc1 = _wf(remove_fourfour=1, add441_amp=AMP, add441_phi=PHI)
        qnm_par = get_qnm_freqtau(['441'], **COMMON)
        conv = COMMON['distance'] * 1e6 * lal.PC_SI / \
            (qnm_par.final_mass * lal.MRSUN_SI)
        n = len(hp0)
        i0 = int(np.floor(float(TOFFSET - hp0.start_time) * hp0.sample_rate))
        # exact slice times: reuse the series' own grid; the generator's
        # time_slice is end-exclusive, so the last sample is not touched
        t = hp0.sample_times.numpy()[i0:-1]
        s = -(AMP / conv) * np.exp(1j * PHI) \
            * np.exp(-1j * qnm_par.omega('441') * (t - TOFFSET))
        exp_p = np.zeros(n)
        exp_c = np.zeros(n)
        exp_p[i0:-1], exp_c[i0:-1] = _expected_pol44(s)
        scale = np.abs(s).max()
        self.assertTrue(np.allclose(hp1.numpy() - hp0.numpy(), exp_p,
                                    atol=1e-8 * scale))
        self.assertTrue(np.allclose(hc1.numpy() - hc0.numpy(), exp_c,
                                    atol=1e-8 * scale))

    def test_add441_linearity_and_phase(self):
        hp0, _ = _wf(remove_fourfour=1)
        d1 = _wf(remove_fourfour=1, add441_amp=AMP)[0].numpy() - hp0.numpy()
        d2 = _wf(remove_fourfour=1, add441_amp=2 * AMP)[0].numpy() - hp0.numpy()
        dpi = _wf(remove_fourfour=1, add441_amp=AMP,
                  add441_phi=np.pi)[0].numpy() - hp0.numpy()
        scale = np.abs(d1).max()
        self.assertTrue(np.allclose(d2, 2 * d1, atol=1e-10 * scale))
        self.assertTrue(np.allclose(dpi, -d1, atol=1e-8 * scale))

    def test_zero_amp_matches_nofourfour(self):
        hp0, _ = _wf(remove_fourfour=1)
        hp1, _ = _wf(remove_fourfour=1, add441_amp=0.0, add441_phi=1.3)
        self.assertTrue(np.array_equal(hp0.numpy(), hp1.numpy()))

    def test_add441_on_removeq_also_works(self):
        # the free 441 composes with the ordinary removeq model too
        hp0, _ = _wf()
        hp1, _ = _wf(add441_amp=AMP)
        self.assertFalse(np.array_equal(hp0.numpy(), hp1.numpy()))


if __name__ == '__main__':
    unittest.main()
