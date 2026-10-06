import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools/abaqus'))
from embedded_host_map import inverse_hex


class HostInterpolationTest(unittest.TestCase):
    def test_rotated_skew_hex_reproduces_affine_displacement(self):
        # Non-diagonal physical Jacobian exposes an incorrect transpose in
        # inverse mapping; use an independent affine reference transformation.
        reference = np.array(((-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
                              (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)))
        mapping = np.array(((0.3, -0.1, 0.04), (0.2, 0.4, 0.05), (-0.05, 0.03, 0.2)))
        offset = np.array((0.02, -0.8, 6.3))
        corners = reference @ mapping.T + offset
        point = mapping @ np.array((0.25, -0.35, 0.6)) + offset
        weights = inverse_hex(point, corners)
        self.assertIsNotNone(weights)
        np.testing.assert_allclose(weights @ corners, point, atol=1e-12)
        gradient = np.array(((2, 3, -1), (0.2, -0.5, 4), (-0.7, 2, 1)))
        values = corners @ gradient.T + np.array((1, 2, 3))
        np.testing.assert_allclose(weights @ values, gradient @ point + (1, 2, 3), atol=1e-12)
        self.assertIsNone(inverse_hex(mapping @ np.array((1.01, 0, 0)) + offset, corners))


if __name__ == '__main__':
    unittest.main()
