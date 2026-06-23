import math
import unittest

import astropy.units as u
import numpy as np
from astropy.cosmology import FlatLambdaCDM, z_at_value
from scipy import constants
from scipy.integrate import quad

from pytgr.massive_graviton import effective_distance, mg_phase_correction, mg_to_lambda_g


class MassiveGravitonFormulaTests(unittest.TestCase):
    def test_effective_distance_matches_will_eds_limit(self):
        cosmo = FlatLambdaCDM(H0=67.4, Om0=1.0, Tcmb0=0.0)
        c_km_s = constants.c / 1000.0

        for z in [0.01, 0.1, 1.0]:
            luminosity_distance = cosmo.luminosity_distance(z).value
            actual = effective_distance(luminosity_distance, cosmology=cosmo)
            expected = (
                2.0
                * c_km_s
                * (1.0 + z)
                * (1.0 - (1.0 + z) ** (-2.5))
                / (5.0 * cosmo.H0.value)
            )
            self.assertAlmostEqual(actual / expected, 1.0, places=6)

    def test_effective_distance_uses_one_cosmology_consistently(self):
        cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=0.0)
        z = 0.25
        luminosity_distance = cosmo.luminosity_distance(z).value

        actual = effective_distance(luminosity_distance, cosmology=cosmo)
        expected = (1.0 + z) * constants.c / 1000.0 / cosmo.H0.value
        expected *= _integrate_will_distance_factor(z, cosmo)
        self.assertAlmostEqual(actual / expected, 1.0, places=7)

    def test_phase_correction_matches_will_absolute_one_over_f_term(self):
        cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=0.0)
        distance_mpc = 1000.0
        frequencies = np.array([20.0, 100.0, 512.0])
        mg = 1e-22

        z = z_at_value(cosmo.luminosity_distance, distance_mpc * u.Mpc).value
        distance_m = effective_distance(distance_mpc, cosmology=cosmo) * 1e6 * constants.parsec
        lambda_g_m = mg_to_lambda_g(mg) * 1000.0
        coefficient = math.pi * distance_m * constants.c / lambda_g_m**2 / (1.0 + z)

        actual = mg_phase_correction(mg, distance_mpc, frequencies, cosmology=cosmo)
        expected = -coefficient / frequencies
        np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=0.0)

    def test_will_beta_maps_to_ppe_coefficient(self):
        chirp_mass_seconds = 35.0 * 4.925490947e-6
        ppe_coefficient = 123.4

        beta = math.pi * chirp_mass_seconds * ppe_coefficient
        frequency = 100.0
        will_phase = -beta / (math.pi * chirp_mass_seconds * frequency)
        ppe_phase = -ppe_coefficient / frequency

        self.assertAlmostEqual(will_phase, ppe_phase, places=15)


def _integrate_will_distance_factor(z, cosmo):
    value, _ = quad(lambda zp: 1.0 / ((1.0 + zp) ** 2 * cosmo.efunc(zp)), 0.0, z)
    return value


if __name__ == "__main__":
    unittest.main()
