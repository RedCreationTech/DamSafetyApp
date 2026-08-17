import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[2]
TOOL_PATH = PROJECT_DIR / 'tools' / 'abaqus' / 'gen_dam_two_step_case.py'
SPEC = importlib.util.spec_from_file_location('gen_dam_two_step_case', TOOL_PATH)
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)


class GenerateDamTwoStepCaseTest(unittest.TestCase):
    def test_generates_static_initialization_and_full_dynamic_step(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            for name in ('mesh.e', 'mass_x.csv', 'mass_y.csv'):
                (work / name).write_text('placeholder\n', encoding='utf-8')
            report = {
                'blocks': {'DAM': {'material': 'CONCRETE'}},
                'nodesets': {'BASE': 2, 'ACCEL': 2, 'POINT_MASS': 1},
                'sidesets': {'WATER': 1},
                'materials': {'concrete': {
                    'elastic': [[30e9, 0.2]], 'density': [[2260]],
                    'damping': {'alpha': 1.95, 'beta': 0.00113}}},
                'amplitudes': {'EQ': [[0.0, -1.0], [2.0, 1.0]]},
                'amplitude_options': {'EQ': {
                    'name': 'EQ', 'time': 'TOTAL TIME'}},
                'initial_boundaries': [{
                    'set': 'BASE', 'dof1': 1, 'dof2': 2, 'value': 0.0,
                    'amplitude': None}],
                'steps': [
                    {'name': 'GRA', 'static': [0.1, 1.0, 1e-15, 1.0],
                     'loads': [
                         {'region': '', 'type': 'GRAV', 'value': 9.8,
                          'parameters': [0.0, -1.0]},
                         {'surface': 'WATER', 'type': 'HP',
                          'value': 592116.0, 'parameters': [60.42, 0.0]}]},
                    {'name': 'EQ', 'dynamic': [0.01, 50.0, 5e-7],
                     'boundaries': [
                         {'set': 'ACCEL', 'dof1': 1, 'dof2': 1,
                          'value': 3.44, 'amplitude': 'EQ',
                          'type': 'ACCELERATION'},
                         {'set': 'BASE', 'dof1': 2, 'dof2': 2,
                          'value': 0.0, 'amplitude': None}]},
                ],
            }
            report_path = work / 'report.json'
            report_path.write_text(json.dumps(report), encoding='utf-8')
            output = work / 'generated'

            static_path, dynamic_path, acceleration_path = GENERATOR.generate_case(
                report_path, work / 'mesh.e', work / 'mass_x.csv',
                work / 'mass_y.csv', output)

            static_text = static_path.read_text(encoding='utf-8')
            dynamic_text = dynamic_path.read_text(encoding='utf-8')
            self.assertIn('type = Steady', static_text)
            self.assertIn('initial_from_file_var = disp_x', dynamic_text)
            self.assertIn('mass_damping_coefficient = 1.95', dynamic_text)
            self.assertIn('stiffness_damping_coefficient = 0.00113',
                          dynamic_text)
            self.assertIn('start_time = 1', dynamic_text)
            self.assertIn('end_time = 51', dynamic_text)
            self.assertIn('time_step_interval = 10', dynamic_text)
            self.assertEqual(acceleration_path.read_text(encoding='utf-8'),
                             '0,-3.44\n2,3.44\n')


if __name__ == '__main__':
    unittest.main()
