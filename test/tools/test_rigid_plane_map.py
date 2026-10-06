import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools/abaqus'))
import abaqus2exodus as mesh
from rigid_plane_map import surface_integrals


class RigidPlaneIntegrationTest(unittest.TestCase):
    def test_exact_moment_for_affine_pressure_on_offcentre_plane(self):
        # The face is deliberately off-centre with respect to the RP. Check
        # integrated forces/moments against analytical polynomial integrals.
        part = mesh.Part('CUBE')
        part.elems = {1: list(range(1, 9))}
        part.elem_types = {1: 'C3D8R'}
        model = mesh.InpModel()
        model.parts = {'CUBE': part}
        model.instances = [mesh.Instance('A', 'CUBE')]
        model.asm_elsets = {'FACE': {'A': [1]}}
        model.surfaces = {'TOP': [('FACE', 'S2')]}
        gm = mesh.GlobalMesh(1e-9, merge_coincident=False)
        coordinates = ((-2, -1, 0), (2, -1, 0), (2, 1, 0), (-2, 1, 0),
                       (-2, -1, 3), (2, -1, 3), (2, 1, 3), (-2, 1, 3))
        for label, point in enumerate(coordinates, 1):
            gm.add(('A', label), point)
        reference = np.array((0.25, -0.5, 3))
        weights, count = surface_integrals(model, gm, 'TOP', reference)
        self.assertEqual(count, 1)
        self.assertEqual(len(weights), 4)
        np.testing.assert_allclose(sum(weights.values()), (8, -2, 4, 0), atol=1e-12)
        # p = 10 + 2*x + 3*y; integrals x*p=2*A*4^2/12,
        # y*p=3*A*2^2/12, then subtract RP offset times resultant.
        total_force = 0.0
        first = np.zeros(3)
        for node, row in weights.items():
            x, y, _ = gm.coords[node - 1]
            pressure = 10 + 2*x + 3*y
            total_force += pressure * row[0]
            first += pressure * row[1:]
        self.assertAlmostEqual(total_force, 80)
        np.testing.assert_allclose(first, (2*8*16/12 - 0.25*80,
                                          3*8*4/12 + 0.5*80, 0), atol=1e-12)
        # Work conjugacy: integral lambda.(theta x r) =
        # theta.integral(r x lambda), using the same nodal moment integrals.
        theta = np.array((0.03, -0.02, 0.04))
        traction = np.array((5.0, -10.0, 3.0))
        integral_r = sum((row[1:] for row in weights.values()), np.zeros(3))
        work = sum(np.dot(traction, np.cross(theta, row[1:])) for row in weights.values())
        self.assertAlmostEqual(work, np.dot(theta, np.cross(integral_r, traction)))


if __name__ == '__main__':
    unittest.main()
