import pathlib
import sys
import types
import unittest
from unittest.mock import patch


class _WaveformStub:
    sample_frequencies = [0.0, 20.0]

    def __getitem__(self, key):
        return self

    def __setitem__(self, key, value):
        pass

    def __imul__(self, other):
        return self


class BaseGRApproximantTests(unittest.TestCase):
    def test_ppe_uses_base_gr_approximant_for_gr_backend(self):
        calls = []
        pycbc = types.ModuleType("pycbc")
        waveform = types.ModuleType("pycbc.waveform")

        def get_fd_waveform(**kwargs):
            calls.append(kwargs.copy())
            return _WaveformStub(), _WaveformStub()

        waveform.get_fd_waveform = get_fd_waveform

        with patch.dict(sys.modules, {"pycbc": pycbc, "pycbc.waveform": waveform}):
            from pytgr import ppe

            hp, hc = ppe.gen_ppe_waveform(
                approximant="ppe",
                base_gr_approximant="SEOBNRv5_ROM",
                mass1=10,
                mass2=10,
                ppebeta2=1e-3,
            )

        self.assertIsInstance(hp, _WaveformStub)
        self.assertIsInstance(hc, _WaveformStub)
        self.assertEqual(calls[0]["approximant"], "SEOBNRv5_ROM")
        self.assertNotIn("base_gr_approximant", calls[0])
        self.assertNotIn("ppebeta2", calls[0])
        self.assertEqual(calls[0]["mass1"], 10)
        self.assertEqual(calls[0]["mass2"], 10)

    def test_ppe_requires_base_gr_approximant(self):
        pycbc = types.ModuleType("pycbc")
        waveform = types.ModuleType("pycbc.waveform")
        waveform.get_fd_waveform = lambda **kwargs: (_WaveformStub(), _WaveformStub())

        with patch.dict(sys.modules, {"pycbc": pycbc, "pycbc.waveform": waveform}):
            from pytgr import ppe

            with self.assertRaisesRegex(ValueError, "base_gr_approximant"):
                ppe.gen_ppe_waveform(approximant="ppe")

    def test_length_hook_uses_base_gr_approximant(self):
        calls = []
        pycbc = types.ModuleType("pycbc")
        waveform_package = types.ModuleType("pycbc.waveform")
        waveform_module = types.ModuleType("pycbc.waveform.waveform")

        def get_waveform_filter_length_in_time(**kwargs):
            calls.append(kwargs.copy())
            return 1.25

        waveform_module.get_waveform_filter_length_in_time = get_waveform_filter_length_in_time

        with patch.dict(
            sys.modules,
            {
                "pycbc": pycbc,
                "pycbc.waveform": waveform_package,
                "pycbc.waveform.waveform": waveform_module,
            },
        ):
            import pytgr

            length = pytgr.length_in_time(
                approximant="ppe",
                base_gr_approximant="SEOBNRv5_ROM",
                mass1=30,
            )

        self.assertEqual(length, 1.25)
        self.assertEqual(calls[0]["approximant"], "SEOBNRv5_ROM")
        self.assertEqual(calls[0]["mass1"], 30)
        self.assertNotIn("base_gr_approximant", calls[0])

    def test_plugin_sources_do_not_reference_legacy_base_key(self):
        package_dir = pathlib.Path(__file__).resolve().parents[1] / "pytgr"
        for path in package_dir.glob("*.py"):
            self.assertNotIn("baseapprox", path.read_text())


if __name__ == "__main__":
    unittest.main()
