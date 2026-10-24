#include "AqueductHHTScheme.h"
#include "gtest/gtest.h"
#include <array>

TEST(AqueductHHTScheme, NormalizationPreservesPhysicalHHTEquation)
{
  for(double alpha:{0.,-.05,-.333333333333,-.41421,-.5})
    for(double mass:{.2,3.,1000.})
    {
      const double a=.7,force=4.2,old=-3.1;
      const auto scaled=AqueductHHTScheme::normalizedResidual(mass*a,force,old,alpha);
      EXPECT_NEAR((1+alpha)*scaled,mass*a+(1+alpha)*force-alpha*old,1e-12);
      EXPECT_GE(AqueductHHTScheme::beta(alpha),.25);
      EXPECT_GE(AqueductHHTScheme::gamma(alpha),.5);
    }
}
TEST(AqueductHHTScheme, AlphaZeroAndConsistentTangent)
{
  EXPECT_DOUBLE_EQ(AqueductHHTScheme::massScale(0),1.);
  const double alpha=-.41421,dt=.007,beta=AqueductHHTScheme::beta(alpha);
  const double M=2.3,K=19.,old_force=7.,predictor=.02,u=.04,h=1e-7;
  auto residual=[&](double x){return AqueductHHTScheme::normalizedResidual(M*(x-predictor)/(beta*dt*dt),K*x,old_force,alpha);};
  EXPECT_NEAR((residual(u+h)-residual(u-h))/(2*h),M/(beta*dt*dt*(1+alpha))+K,1e-5);
  EXPECT_THROW(AqueductHHTScheme::massScale(-1),std::invalid_argument);
  EXPECT_THROW(AqueductHHTScheme::massScale(.1),std::invalid_argument);
}
