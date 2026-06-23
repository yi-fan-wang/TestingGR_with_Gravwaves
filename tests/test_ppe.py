import math
import os
import unittest

import numpy as np
from pycbc import conversions

from pytgr import ppe


os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-pytgr-tests")


class PpePhaseTests(unittest.TestCase):
    def test_ppe_beta_exponent_matches_pn_index(self):
        self.assertAlmostEqual(ppe.ppe_beta_exponent(1), -4.0 / 3.0)
        self.assertAlmostEqual(ppe.ppe_beta_exponent(2), -1.0)
        self.assertAlmostEqual(ppe.ppe_beta_exponent(7), 2.0 / 3.0)

    def test_phase_shift_uses_original_ppe_beta_terms(self):
        frequencies = np.array([20.0, 50.0, 100.0])
        chirp_mass_seconds = 35.0 * ppe.MTSUN_SI

        actual = ppe._ppe_phase_shift(
            frequencies,
            ppe_beta_terms=[(2, 0.25)],
            chirp_mass_seconds=chirp_mass_seconds,
        )
        expected = 0.25 / (np.pi * chirp_mass_seconds * frequencies)

        np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0.0)

    def test_component_masses_define_detector_chirp_mass(self):
        actual = conversions.mchirp_from_mass1_mass2(10.0, 10.0)
        expected = 10.0 / 2.0**(1.0 / 5.0)
        self.assertAlmostEqual(actual, expected)

    def test_ppebeta2_waveform_ratio_matches_beta_over_u(self):
        from pycbc.waveform import get_fd_waveform

        beta2 = 1e-3
        kwargs = dict(
            approximant="ppe",
            base_gr_approximant="TaylorF2",
            mass1=10,
            mass2=10,
            spin1z=0,
            spin2z=0,
            distance=500,
            inclination=0,
            delta_f=1.0,
            f_lower=20.0,
            ppebeta2=beta2,
        )
        hp_ppe, _ = ppe.gen_ppe_waveform(**kwargs.copy())

        gr_kwargs = kwargs.copy()
        gr_kwargs.pop("approximant")
        gr_kwargs.pop("base_gr_approximant")
        gr_kwargs.pop("ppebeta2")
        hp_gr, _ = get_fd_waveform(approximant="TaylorF2", **gr_kwargs)

        chirp_mass_seconds = conversions.mchirp_from_mass1_mass2(
            kwargs["mass1"], kwargs["mass2"]
        ) * ppe.MTSUN_SI
        for frequency in [20, 50, 100, 200]:
            idx = int(frequency / hp_ppe.delta_f)
            ratio = complex(hp_ppe[idx] / hp_gr[idx])
            expected_phase = beta2 / (np.pi * chirp_mass_seconds * frequency)
            actual_phase = math.atan2(ratio.imag, ratio.real)
            self.assertAlmostEqual(actual_phase, expected_phase, places=12)


if __name__ == "__main__":
    unittest.main()
