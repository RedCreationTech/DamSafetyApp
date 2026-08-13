import importlib.util
import unittest
from pathlib import Path

import numpy as np


PROJECT_DIR = Path(__file__).resolve().parents[2]
TOOL_PATH = PROJECT_DIR / 'tools' / 'abaqus' / 'beam_field_solver.py'
SPEC = importlib.util.spec_from_file_location('beam_field_solver', TOOL_PATH)
SOLVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SOLVER)


class BeamFieldSolverTest(unittest.TestCase):
    def test_circular_section_limits(self):
        self.assertEqual(
            SOLVER.section_fiber_limits('CIRC', [25.0]),
            (25.0, 25.0))

    def test_axial_stress_recovery(self):
        elastic_modulus = 200.0
        area = 2.0
        length = 1.0
        iy = iz = j = 1.0
        stiffness = SOLVER.beam_ke_local(
            elastic_modulus, 80.0, area, iy, iz, j, length, area, area)
        displacement = np.zeros((1, 2, 6))
        displacement[0, 1, 0] = 0.01
        records = {'beam': [{
            'nodes': (0, 1),
            'transform': np.eye(12),
            'ke': stiffness,
            'keep': list(range(12)),
            'area': area,
            'iy': iy,
            'iz': iz,
            'j': j,
            'cy': 1.0,
            'cz': 1.0,
        }]}

        fields = SOLVER.recover_beam_stress(displacement, records)

        self.assertAlmostEqual(fields['axial_stress']['beam'][0, 0], 2.0)
        self.assertAlmostEqual(fields['bending_stress']['beam'][0, 0], 0.0)
        self.assertAlmostEqual(fields['torsional_shear']['beam'][0, 0], 0.0)
        self.assertAlmostEqual(fields['vonmises_stress']['beam'][0, 0], 2.0)


if __name__ == '__main__':
    unittest.main()
