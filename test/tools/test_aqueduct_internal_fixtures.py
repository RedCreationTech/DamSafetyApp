import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/abaqus'))
from build_aqueduct_internal_fixtures import hex_matrices,steel_oracle,XYZ

class InternalPhysicsTest(unittest.TestCase):
    def test_independent_hex8_mass_and_affine_energy(self):
        M,K=hex_matrices();xyz=np.array(XYZ);E=30e9;nu=.2
        np.testing.assert_allclose(M,M.T,atol=1e-14);np.testing.assert_allclose(K,K.T,atol=2e-6)
        for j in range(3):
            u=np.zeros((8,3));u[:,j]=1;np.testing.assert_allclose(K@u.ravel(),0,atol=1e-5)
            self.assertAlmostEqual(float(u.ravel()@M@u.ravel()),2640,places=10)
        u=xyz*np.array([.001,-nu*.001,-nu*.001]);self.assertAlmostEqual(float(u.ravel()@K@u.ravel()),E*.001**2,places=7)
        vals=np.linalg.eigvalsh(K);self.assertEqual(int(np.sum(vals>1e-3)),18)
    def test_source_steel_end_plateau_and_reverse_hardening(self):
        r=steel_oracle([0,.002,.12,.12-696e6/200e9,-.12,0])
        self.assertAlmostEqual(r[1]['sxx'],400e6,delta=1e-5)
        self.assertAlmostEqual(r[2]['sxx'],696e6,delta=1e-4);self.assertGreater(r[2]['peq'],.1)
        self.assertAlmostEqual(r[3]['sxx'],0,delta=1e-4)
        self.assertAlmostEqual(r[4]['sxx'],-696e6,delta=1e-4)
        self.assertGreater(r[4]['peq'],r[2]['peq'])
if __name__=='__main__':unittest.main()
