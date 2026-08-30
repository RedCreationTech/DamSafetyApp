import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools' / 'abaqus'))
from compare_uniaxial_tension_diagnostic import metrics


class ComparisonMetricsTest(unittest.TestCase):
    def test_zero_reference_retains_absolute_error_without_division(self):
        result = metrics(np.zeros(3), np.array([0., .1, 0.]), 0.)
        self.assertEqual(result['max_abs_error'], .1)
        self.assertEqual(result['worst_element'], 2)
        self.assertIsNone(result['max_local_relative_percent'])
        self.assertIsNone(result['nrmse_global_peak_percent'])

    def test_global_normalization_does_not_mask_pointwise_failure(self):
        result = metrics(np.array([.01, 1.]), np.array([.02, 1.]), 1.)
        self.assertLess(result['nrmse_global_peak_percent'], 5)
        self.assertEqual(result['max_local_relative_percent'], 100.)
        self.assertEqual(result['within_5_percent_of_local_reference_fraction'], .5)

    def test_near_zero_points_excluded_only_from_relative_metric(self):
        result = metrics(np.array([0., 1e-5, 1.]), np.array([.02, .02, 1.]), 1.)
        self.assertEqual(result['relative_eligible_elements'], 1)
        self.assertEqual(result['max_abs_error'], .02)
        self.assertEqual(result['within_5_percent_of_local_reference_fraction'], 1.)


if __name__ == '__main__':unittest.main()
