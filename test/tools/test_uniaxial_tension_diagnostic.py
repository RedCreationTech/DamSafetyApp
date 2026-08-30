import sys
import tempfile
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / 'tools' / 'abaqus'
sys.path.insert(0, str(TOOLS))
import build_uniaxial_tension_diagnostic as builder

MATERIAL = '''*Material, name=concrete
*Elastic
30000, 0.2
*Concrete Damaged Plasticity
36, 0.1, 1.16, 0.667, 0.0005
*Concrete Compression Hardening
20, 0
22, 0.001
*Concrete Compression Damage
0, 0
0.1, 0.001
*Concrete Tension Stiffening
3, 0
2, 0.001
*Concrete Tension Damage
0, 0
0.1, 0.001
'''


class TensionDiagnosticTest(unittest.TestCase):
    def test_mpa_conversion_does_not_change_strain_damage_or_time(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'material.inp'
            source.write_text(MATERIAL)
            raw, si = builder.material_si(source)
            self.assertEqual(si['elastic']['youngs_modulus_pa'], 3e10)
            self.assertEqual(si['tables']['tension_stiffening'][0], (3e6, 0))
            self.assertEqual(si['tables']['compression_hardening'][1], (22e6, .001))
            self.assertEqual(si['tables']['tension_damage'], raw['tables']['tension_damage'])
            self.assertEqual(si['plasticity'], raw['plasticity'])
            self.assertEqual(si['recovery'], {'tension_recovery': 0, 'compression_recovery': 1})
            self.assertEqual(source.read_text(), MATERIAL)
            for mode in ('tension', 'compression'):
                for a, b in zip(raw['derived'][mode], si['derived'][mode]):
                    self.assertAlmostEqual(a['equivalent_plastic_strain'], b['equivalent_plastic_strain'])

    def test_source_mismatch_rejected_before_output_creation(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'changed.inp'
            source.write_text(MATERIAL)
            output = Path(tmp) / 'out'
            with self.assertRaisesRegex(ValueError, 'Source changed'):
                builder.build(source, output)
            self.assertFalse(output.exists())

    def test_loading_is_z_ramp_and_bottom_is_transversely_free(self):
        text = builder.moose_input({'elastic': {'youngs_modulus_pa': 3e10, 'poissons_ratio': .2}},
                                   '[Materials]\n  [cdp_stress_update]\n  []\n[]\n')
        self.assertIn("expression = '2.5e-5*t'", text)
        self.assertIn('dt = 0.01', text)
        self.assertIn('end_time = 1', text)
        self.assertIn('dtmin = 1e-15', text)
        self.assertNotIn('[bottom_x]', text)
        self.assertNotIn('[bottom_y]', text)
        self.assertIn('[top_x_gauge]', text)
        self.assertIn('[top_y_gauge]', text)


if __name__ == '__main__':
    unittest.main()
