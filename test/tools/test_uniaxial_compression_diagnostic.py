import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools' / 'abaqus'))
import build_uniaxial_compression_diagnostic as builder


class CompressionDiagnosticTest(unittest.TestCase):
    def test_compression_ramp_sign_magnitude_and_output(self):
        text = builder.moose_input({'elastic': {'youngs_modulus_pa': 2.97915e10, 'poissons_ratio': .2}},
                                   '[Materials]\n  [cdp_stress_update]\n  []\n[]\n')
        for value in ("expression = '-0.0025*t'", 'end_time = 1', 'dt = 0.01', 'dtmin = 1e-15',
                      'nl_rel_tol = 1e-9', 'nl_abs_tol = 1e-8', '[min_stress_zz]', 'value_type = min',
                      'file_base = uniaxial_compression', builder.SOURCE_SHA):
            self.assertIn(value, text)
        for value in ('C3D8R', 'uniaxial_tension', '[bottom_x]', '[bottom_y]', '2.5e-5*t'):
            self.assertNotIn(value, text)

    def test_changed_source_rejected_without_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'compression.inp'; source.write_text('changed source')
            output = Path(tmp) / 'out'
            with self.assertRaisesRegex(ValueError, 'Source changed'):
                builder.build(source, output)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
