import unittest
import numpy as np
from aqueduct_s2s import mortar

class MortarTest(unittest.TestCase):
    def test_nonmatching_affine_and_virtual_work(self):
        # A coarse secondary surface overlaps four refined primary quads.
        p=np.array([[x,y,0.] for y in [0.,.5,1.] for x in [0.,.5,1.]]+[[0,0,0],[1,0,0],[1,1,0],[0,1,0]])
        secondary=[dict(nodes=[10,11,12,13])]
        primary=[dict(nodes=[a,a+1,a+4,a+3]) for a in (1,2,4,5)]
        rows,meta=mortar(p,secondary,primary,[10,11,12,13])
        self.assertAlmostEqual(meta['coverage_fraction'],1)
        C=np.zeros((4,len(p)))
        for i,row in enumerate(rows):
            for t in row['scalar_terms']:C[i,t['node']]=t['coefficient']
        self.assertLess(np.max(abs(C@np.column_stack([np.ones(len(p)),p[:,:2]]))),1e-12)
        self.assertTrue(all(len(r['primary_nodes'])==9 for r in rows))
        lm=np.array([1.,2.,-3.,4.]);force=C.T@lm
        self.assertLess(abs(sum(force)),1e-12)
        self.assertLess(np.linalg.norm(p.T@force),1e-12)
        du=np.random.default_rng(12).normal(size=len(p))
        self.assertAlmostEqual(force@du,lm@(C@du),places=12)
        rows2,_=mortar(p,secondary,primary,[10,11,12,13],order=7)
        for a,b in zip(rows,rows2):np.testing.assert_allclose([t['coefficient'] for t in a['scalar_terms']],[t['coefficient'] for t in b['scalar_terms']],atol=1e-13)
    def test_unbound_nodes_not_forced(self):
        p=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,0],[1,0,0],[1,1,0],[0,1,0]])
        with self.assertRaisesRegex(ValueError,'fully bound'):mortar(p,[dict(nodes=[5,6,7,8])],[dict(nodes=[1,2,3,4])],[5,6,7])

if __name__=='__main__':unittest.main()
