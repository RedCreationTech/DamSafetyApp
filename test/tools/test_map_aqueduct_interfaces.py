import sys
from pathlib import Path
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/abaqus'))
import map_aqueduct_interfaces as tool

class InterfaceMappingTest(unittest.TestCase):
    def setUp(self):self.q=np.array([[0.,0,0],[1,0,0],[1.3,1,0],[.3,1,0]])
    def test_bounded_projection_and_free_edge_stationarity(self):
        target=np.array([1.5,.4,.01]);p=tool.project_quad(target,self.q)
        self.assertTrue((p['weights']>=-1e-15).all())
        self.assertAlmostEqual(p['weights'].sum(),1)
        self.assertAlmostEqual(p['natural'][0],1)
        self.assertAlmostEqual(float((target-p['point'])@(self.q[2]-self.q[1])),0,places=12)
    def test_rotation_and_offset_translation_reproduce_rigid_motion(self):
        p=tool.project_quad(np.array([.6,.4,.2]),self.q)
        A=tool.rotation_operator(p['gradients'],p['normal']);r=np.array([0,0,.2])
        P=np.array([N*np.eye(3)-tool.cross_matrix(r)@a for N,a in zip(p['weights'],A)])
        theta=np.array([.02,-.03,.04]);translation=np.array([.1,.2,-.1])
        u=np.cross(theta,self.q)+translation
        np.testing.assert_allclose(np.einsum('ijk,ik->j',A,u),theta,atol=1e-14)
        np.testing.assert_allclose(np.einsum('ijk,ik->j',P,u),np.cross(theta,p['point']+r)+translation,atol=1e-14)
    def test_conjugate_transfer_balances_force_and_moment(self):
        p=tool.project_quad(np.array([.6,.4,.2]),self.q);r=np.array([0,0,.2]);A=tool.rotation_operator(p['gradients'],p['normal'])
        P=np.array([N*np.eye(3)-tool.cross_matrix(r)@a for N,a in zip(p['weights'],A)])
        force=np.array([12.,-7,4]);moment=np.array([3.,5,-2])
        host=-np.einsum('ijk,j->ik',P,force)-np.einsum('ijk,j->ik',A,moment)
        np.testing.assert_allclose(host.sum(0)+force,0,atol=1e-12)
        np.testing.assert_allclose(np.cross(self.q,host).sum(0)+np.cross(p['point']+r,force)+moment,0,atol=1e-12)
    def test_projection_does_not_fall_back_to_unrelated_surface(self):
        facets=[dict(nodes=[1,2,3,4],instance='HOST',element_label=1,face='S2')]
        with self.assertRaises(ValueError):tool.project_surface([np.array([.5,.5,.06])],facets,self.q)

if __name__=='__main__':unittest.main()
