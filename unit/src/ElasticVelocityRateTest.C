#include "ElasticVelocityRate.h"
#include "gtest/gtest.h"
#include <cmath>
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

TEST(ElasticVelocityRate, ConstantAccelerationUsesAcceptedHistoryOnRetry)
{
  const double u0 = .02, v0 = .07, a0 = -.3, beta = .275625, gamma = .55;
  for (double h : {.01, .005, .0037})
  {
    const double u = u0 + h * v0 + .5 * h * h * a0;
    EXPECT_NEAR(ElasticVelocityRate::trialVelocity(u,u0,v0,a0,h,beta,gamma), v0 + h * a0, 1e-13);
  }
  // A discarded trial cannot become old_u/old_v/old_a for a smaller retry.
  const double rejected = ElasticVelocityRate::trialVelocity(.031,u0,v0,a0,.01,beta,gamma);
  const double retry_u = u0 + .005 * v0 + .5 * .005 * .005 * a0;
  EXPECT_GT(std::abs(ElasticVelocityRate::trialVelocity(retry_u,.031,rejected,a0,.005,beta,gamma) - (v0 + .005*a0)), 1);
}

TEST(ElasticVelocityRate, MixedElasticJacobianHasDistinctThicknessColumns)
{
  const double beta = .275625, gamma = .55, alpha = -.05, zeta = .03;
  const double kuu = 11, kuz = 3, kzu = 3, kzz = 7;
  for (double h : {.01, .005, .0037})
  {
    auto force = [&](double u, double z, double r) {
      const double v = ElasticVelocityRate::trialVelocity(u,.01,.07,-.3,h,beta,gamma);
      return (1+alpha)*(kuu*u+kuz*z)-alpha*(kuu*.01+kuz*.02) +
          zeta*((1+alpha)*(kuu*v+kuz*r)-alpha*(kuu*.07+kuz*.03));
    };
    auto constraint = [&](double u, double r) {
      return kzu*ElasticVelocityRate::trialVelocity(u,.01,.07,-.3,h,beta,gamma)+kzz*r;
    };
    const double e = 1e-6, c = ElasticVelocityRate::displacementDerivative(h,beta,gamma);
    EXPECT_NEAR((force(.013+e,.02,.03)-force(.013-e,.02,.03))/(2*e), (1+alpha)*kuu*(1+zeta*c), 1e-7);
    EXPECT_NEAR((force(.013,.02+e,.03)-force(.013,.02-e,.03))/(2*e), (1+alpha)*kuz, 1e-7);
    EXPECT_NEAR((force(.013,.02,.03+e)-force(.013,.02,.03-e))/(2*e), zeta*(1+alpha)*kuz, 1e-7);
    EXPECT_NEAR((constraint(.013+e,.03)-constraint(.013-e,.03))/(2*e), kzu*c, 1e-7);
    EXPECT_NEAR((constraint(.013,.03+e)-constraint(.013,.03-e))/(2*e), kzz, 1e-7);
    EXPECT_GT(std::abs((1+alpha)*kuu*(1+zeta*c)-(1+alpha)*kuu*(1+zeta/h)), 1);
  }
}

TEST(ElasticVelocityRate, InitialAssemblyDoesNotAttemptZeroIntervalNewmark)
{
  EXPECT_THROW(ElasticVelocityRate::trialVelocity(0.,0.,0.,0.,0.,.25,.5), std::invalid_argument);
  EXPECT_THROW(ElasticVelocityRate::displacementDerivative(.01,0.,.5), std::invalid_argument);
}

TEST(ElasticVelocityRate, ObservableClosureAndSchurDerivativeAgree)
{
  const double E = 100, nu = .2;
  const double mu = E / (2*(1+nu)), lambda = E*nu / ((1+nu)*(1-2*nu));
  RankFourTensor C;
  for (unsigned int i=0;i<3;++i)
    for (unsigned int j=0;j<3;++j)
      for (unsigned int k=0;k<3;++k)
        for (unsigned int l=0;l<3;++l)
          C(i,j,k,l) = lambda*(i==j)*(k==l)+mu*((i==k)*(j==l)+(i==l)*(j==k));
  RankTwoTensor rate;
  rate(0,0)=.07;rate(1,1)=-.02;rate(0,1)=rate(1,0)=.013;
  const auto closed = ElasticVelocityRate::closePlaneStressRate(C,rate);
  const auto tangent = ElasticVelocityRate::planeStressRateTangent(C);
  EXPECT_NEAR((C*closed)(2,2),0,1e-13);
  EXPECT_NEAR((C*closed - tangent*rate).L2norm(),0,1e-13);
  // A plane-strain rate is an explicit wrong alternative, not another valid closure.
  EXPECT_GT(std::abs((C*rate)(0,0)-(C*closed)(0,0)), .1);
  for (double h : {.01,.005,.0037})
  {
    const double c=ElasticVelocityRate::displacementDerivative(h,.275625,.55), e=1e-6;
    RankTwoTensor direction;
    direction(0,0)=.3;direction(1,1)=-.1;
    direction(0,1)=direction(1,0)=.05;
    const auto plus=ElasticVelocityRate::closePlaneStressRate(C,rate+direction*(e*c));
    const auto minus=ElasticVelocityRate::closePlaneStressRate(C,rate-direction*(e*c));
    EXPECT_NEAR(((C*(plus-minus))/(2*e) - (tangent*direction)*c).L2norm(),0,1e-7);
  }
  EXPECT_THROW(ElasticVelocityRate::planeStressRateTangent(RankFourTensor()),std::invalid_argument);
}
