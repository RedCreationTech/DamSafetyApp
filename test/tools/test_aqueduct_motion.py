import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/abaqus'))
from prepare_aqueduct_motion import baseline,integrate

class MotionTest(unittest.TestCase):
    def test_exact_constant_and_linear_acceleration_are_removed(self):
        t=np.array([0.,.013,.2,.79,1.])
        for a in (np.ones(5)*2,2+3*t):
            r=baseline(t,a);np.testing.assert_allclose(r['corrected'],0,atol=2e-10)
    def test_independent_quadrature_stationarity_for_nonuniform_knots(self):
        t=np.array([0.,.03,.2,.5,.77,1.]);a=np.array([.2,.4,-.9,.3,-.1,.6]);r=baseline(t,a)
        v,u=integrate(t,a);nodes,weights=np.polynomial.legendre.leggauss(6);mom=np.zeros(3)
        c=np.array(r['coefficients'])
        for j,h in enumerate(np.diff(t)):
            x=t[j]+(nodes+1)*h/2;z=x-t[j];slope=(a[j+1]-a[j])/h
            vc=v[j]+a[j]*z+slope*z*z/2+c[0]*x+c[1]*x*x/2+c[2]*x*x*x/3
            for k in range(3):mom[k]+=np.sum(weights*vc*x**(k+1)/(k+1))*h/2
        np.testing.assert_allclose(mom,0,atol=1e-12)
    def test_domains_and_raw_inputs_are_not_mutated(self):
        t=np.array([0.,.4,1.]);a=np.array([0.,1.,0.]);old=a.copy();baseline(t,a);np.testing.assert_array_equal(a,old)
        for times,accel in [([0,0],[1,2]),([0,1],[1]),([0,1],[1,float('nan')])]:
            with self.assertRaises(ValueError):baseline(times,accel)
        with self.assertRaises(ValueError):baseline([.1,1],[0,1])
if __name__=='__main__':unittest.main()
