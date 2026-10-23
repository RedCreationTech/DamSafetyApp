import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

TOOLS = Path(__file__).resolve().parents[2] / 'tools' / 'abaqus'
sys.path.insert(0, str(TOOLS))
spec = importlib.util.spec_from_file_location('audit_aqueduct_geometry', TOOLS / 'audit_aqueduct_geometry.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class GeometryAuditTest(unittest.TestCase):
    def setUp(self):
        self.cube = (audit.SIGNS + 1) / 2

    def test_volume_integration_detects_orientation_not_only_connectivity(self):
        mirrored = self.cube.copy(); mirrored[:, 0] *= -1
        det, volume = audit.hex_volumes(np.array([self.cube, mirrored]))
        np.testing.assert_allclose(volume, [1, -1], rtol=1e-14)
        self.assertTrue((det[0] > 0).all())
        self.assertTrue((det[1] < 0).all())

    def test_host_shape_weights_reproduce_an_affine_displacement_field(self):
        point = np.array([0.2, 0.3, 0.4])
        found, unresolved = audit.locate_embedded([(10, point)], self.cube[None], [('A', 7)])
        self.assertFalse(unresolved)
        weights = np.array(found[0]['weights'])
        A = np.array([[2,3,4], [0,4,2], [3,1,2]])
        b = np.array([3,7,2])
        np.testing.assert_allclose(weights @ (self.cube @ A + b), point @ A + b,
                                   atol=1e-13)
        self.assertEqual(found[0]['host'], ['A', 7])

    def test_rotated_skew_host_inverse_uses_physical_by_natural_jacobian(self):
        transform = np.array([[0,-2,0.3], [1,0,0.4], [0.1,0.2,1.5]])
        host = self.cube @ transform.T + [5,7,3]
        natural = np.array([0.2,-0.4,0.6])
        weights = np.array(audit.mesh._hex_shape(*natural))
        point = weights @ host
        recovered = audit.mesh.point_in_hex(point, host, tol=1e-8)
        self.assertIsNotNone(recovered)
        np.testing.assert_allclose(recovered, natural, atol=1e-8)
        _, jac = audit.mesh._hex_map(host, *natural)
        for j in range(3):
            plus = natural.copy(); plus[j] += 1e-6
            minus = natural.copy(); minus[j] -= 1e-6
            p1, _ = audit.mesh._hex_map(host, *plus)
            p2, _ = audit.mesh._hex_map(host, *minus)
            np.testing.assert_allclose(np.array(jac)[:,j],
                                       (np.array(p1)-np.array(p2))/2e-6, atol=1e-8)

    def test_outside_declared_host_is_reported_without_nearest_fallback(self):
        found, missing = audit.locate_embedded([(11, np.array([2.,0.5,0.5]))],
                                               self.cube[None], [('A', 7)])
        self.assertFalse(found)
        self.assertEqual(missing[0]['node'], 11)

    def test_host_boundary_deterministic_selection_uses_source_identity(self):
        hosts = np.array([self.cube, self.cube + [1,0,0]])
        found, missing = audit.locate_embedded([(11,np.array([1.,0.5,0.5]))],
                                               hosts, [('Z',2),('A',1)])
        self.assertFalse(missing)
        self.assertEqual(found[0]['host'], ['A',1])


if __name__ == '__main__':
    unittest.main()
