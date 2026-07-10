import unittest

import numpy as np

try:
    from pycbc.types import TimeSeries
except ImportError:
    TimeSeries = None

try:
    from tgr.nrsurqnm import (QNMTable, least_square_qnmfitting,
                              weighted_least_square_qnmfitting)
except ImportError:
    QNMTable = None
    least_square_qnmfitting = None
    weighted_least_square_qnmfitting = None


def make_three_mode_data():
    '''Synthetic ringdown: three overtones with known complex amplitudes.

    Frequencies/damping times are GW150914-like; the exact values are
    irrelevant, only that the overtones decay progressively faster.
    '''
    qnm_par = QNMTable(final_mass=62.0, final_spin=0.68,
                       freq={'220': 250.7, '221': 248.6, '222': 244.0},
                       tau={'220': 0.00415, '221': 0.00138, '222': 0.00082})
    A_true = {'220': 4.7e-21 * np.exp(0.3j),
              '221': 8.0e-21 * np.exp(-1.2j),
              '222': 9.0e-21 * np.exp(2.0j)}
    delta_t = 1.0 / 16384
    t0 = 0.001
    t = np.arange(600) * delta_t
    data = np.zeros(len(t), dtype=np.complex128)
    for m, A in A_true.items():
        data += A * np.exp(-1j * qnm_par.omega(m) * (t - t0))
    h_target = TimeSeries(data, delta_t=delta_t, epoch=0)
    return qnm_par, A_true, t0, h_target


@unittest.skipUnless(TimeSeries is not None and QNMTable is not None,
                     "pycbc/lal (or tgr optional dependencies) not available")
class WeightedQNMFittingTests(unittest.TestCase):
    def setUp(self):
        self.qnm_par, self.A_true, self.t0, self.h_target = make_three_mode_data()

    def test_unweighted_recovers_complete_model(self):
        # data lies exactly in the span of the three templates
        A_fit, _, _ = least_square_qnmfitting(
            ['220', '221', '222'], self.qnm_par, self.t0, self.h_target)
        for m in self.A_true:
            self.assertLess(abs(A_fit[m] - self.A_true[m]) / abs(self.A_true[m]),
                            1e-3, msg=f"mode {m}")

    def test_uniform_weights_match_unweighted(self):
        fit_modes = ['220', '221']
        A_ols, _, h_slice = least_square_qnmfitting(
            fit_modes, self.qnm_par, self.t0, self.h_target)
        A_uni, _, _ = least_square_qnmfitting(
            fit_modes, self.qnm_par, self.t0, self.h_target,
            weights=np.ones(len(h_slice)))
        for m in fit_modes:
            self.assertLess(abs(A_uni[m] - A_ols[m]) / abs(A_ols[m]), 1e-10)

    def test_oracle_downweighting_beats_ordinary_ls(self):
        # misspecified fit: 222 is in the data but not in the fit model
        fit_modes = ['220', '221']
        A_ols, _, h_slice = least_square_qnmfitting(
            fit_modes, self.qnm_par, self.t0, self.h_target)

        # oracle weights from the true omitted-mode envelope
        t = h_slice.sample_times.numpy()
        sigma2 = np.abs(self.A_true['222'])**2 \
            * np.exp(-2 * (t - self.t0) / self.qnm_par.tau['222'])
        eps = 1e-4 * abs(self.A_true['222'])
        A_wls, _, _ = least_square_qnmfitting(
            fit_modes, self.qnm_par, self.t0, self.h_target,
            weights=1 / (eps**2 + sigma2))

        for m in fit_modes:
            err_ols = abs(A_ols[m] - self.A_true[m])
            err_wls = abs(A_wls[m] - self.A_true[m])
            self.assertLess(err_wls, err_ols, msg=f"mode {m}")

    def test_two_pass_weighted_fit(self):
        fit_modes = ['220', '221']
        A_ols, _, _ = least_square_qnmfitting(
            fit_modes, self.qnm_par, self.t0, self.h_target)
        A_wls, _, h_slice, info = weighted_least_square_qnmfitting(
            fit_modes, self.qnm_par, self.t0, self.h_target,
            omitted_modes=['222'])

        # pass 1 fits the complete model, so the omitted amplitude is accurate
        self.assertLess(abs(info['A_omitted']['222'] - self.A_true['222'])
                        / abs(self.A_true['222']), 1e-3)
        # down-weighting improves the retained-mode amplitudes
        for m in fit_modes:
            err_ols = abs(A_ols[m] - self.A_true[m])
            err_wls = abs(A_wls[m] - self.A_true[m])
            self.assertLess(err_wls, err_ols, msg=f"mode {m}")
        # default floor is 1% of the peak data amplitude in the window
        expected_floor = 1e-2 * np.abs(h_slice.data).max()
        self.assertAlmostEqual(info['epsilon_floor'] / expected_floor, 1.0,
                               places=6)
        self.assertEqual(len(info['weights']), len(h_slice))
        self.assertEqual(len(info['sigma']), len(h_slice))

    def test_explicit_epsilon_floor_is_used(self):
        floor = 5e-23
        _, _, _, info = weighted_least_square_qnmfitting(
            ['220', '221'], self.qnm_par, self.t0, self.h_target,
            omitted_modes=['222'], epsilon_floor=floor)
        self.assertEqual(info['epsilon_floor'], floor)

    def test_invalid_inputs_raise(self):
        with self.assertRaises(ValueError):
            weighted_least_square_qnmfitting(
                ['220', '221'], self.qnm_par, self.t0, self.h_target,
                omitted_modes=['221'])  # overlaps fit_modes
        with self.assertRaises(ValueError):
            weighted_least_square_qnmfitting(
                ['220', '221'], self.qnm_par, self.t0, self.h_target,
                omitted_modes=[])  # nothing to build the envelope from
        with self.assertRaises(ValueError):
            least_square_qnmfitting(
                ['220', '221'], self.qnm_par, self.t0, self.h_target,
                weights=np.ones(3))  # wrong length
        _, _, h_slice = least_square_qnmfitting(
            ['220', '221'], self.qnm_par, self.t0, self.h_target)
        bad = np.ones(len(h_slice))
        bad[0] = -1.0
        with self.assertRaises(ValueError):
            least_square_qnmfitting(
                ['220', '221'], self.qnm_par, self.t0, self.h_target,
                weights=bad)  # non-positive weight


if __name__ == '__main__':
    unittest.main()
