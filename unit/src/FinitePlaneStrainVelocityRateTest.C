#include "FinitePlaneStrainVelocityRate.h"
#include "gtest/gtest.h"
#include "ADUtils.h"
#include <algorithm>

namespace
{
constexpr Real E = 2.97915e10, nu = .2, zeta = .003604;
constexpr Real beta = .275625, gamma = .55, alpha = -.05, dt = .002;
RankTwoTensor tensor(Real xx, Real xy, Real yx, Real yy)
{
  RankTwoTensor a;
  a(0,0)=xx; a(0,1)=xy; a(1,0)=yx; a(1,1)=yy;
  return a;
}
}

TEST(FinitePlaneStrainVelocityRate, ReusesNewmarkAcceptedHistoryAcrossRetries)
{
  const Real u0=.02, v0=.07, a0=-.3;
  for (Real h : {.01, .005, .0037})
  {
    const Real u=u0+h*v0+.5*h*h*a0;
    EXPECT_NEAR(ElasticVelocityRate::trialVelocity(u,u0,v0,a0,h,beta,gamma),v0+h*a0,1e-13);
  }
  const auto before=ElasticVelocityRate::trialVelocity(.031,u0,v0,a0,.002,beta,gamma);
  (void)ElasticVelocityRate::trialVelocity(.035,u0,v0,a0,.004,beta,gamma);
  EXPECT_DOUBLE_EQ(before,ElasticVelocityRate::trialVelocity(.031,u0,v0,a0,.002,beta,gamma));
}

TEST(FinitePlaneStrainVelocityRate, PlaneStrainRetainsNonzeroOutOfPlaneStress)
{
  const auto f=FinitePlaneStrainVelocityRate::deformation(RankTwoTensor());
  const auto rate=tensor(.001,.0002,.0002,-.0003);
  const auto stress=FinitePlaneStrainVelocityRate::stress(rate,f,E,nu);
  const Real lambda=E*nu/((1+nu)*(1-2*nu)), mu=E/(2*(1+nu));
  EXPECT_NEAR(stress(2,2),lambda*.0007,1e-8);
  EXPECT_NEAR(stress(0,0),2*mu*.001+lambda*.0007,1e-8);
  EXPECT_NEAR(stress(0,1),2*mu*.0002,1e-8);
}

TEST(FinitePlaneStrainVelocityRate, InstantRigidRotationHasZeroDampingStress)
{
  for (Real angle : {0.,.3,1.2})
  {
    auto f=tensor(std::cos(angle),-std::sin(angle),std::sin(angle),std::cos(angle));
    f(2,2)=1.;
    const auto omega=tensor(0.,-.7,.7,0.);
    const auto stress=FinitePlaneStrainVelocityRate::stress(omega*f,f,E,nu);
    EXPECT_NEAR(stress.L2norm(),0.,2e-5);
  }
}

TEST(FinitePlaneStrainVelocityRate, StaticPredeformationWithZeroVelocityIsAtRest)
{
  const auto g=tensor(.12,.07,-.02,-.03);
  const auto p=FinitePlaneStrainVelocityRate::nominalStress(
      g,g,RankTwoTensor(),RankTwoTensor(),dt,beta,gamma,alpha,zeta,E,nu);
  EXPECT_DOUBLE_EQ(p.L2norm(),0.);
}

TEST(FinitePlaneStrainVelocityRate, HHTUsesAcceptedOldSpatialRateAndCurrentPullback)
{
  const auto g=tensor(.12,.03,-.02,-.04), old=tensor(.1,.02,-.01,-.03);
  const auto v=tensor(.004,.002,-.003,.001), a=tensor(.01,.02,.03,-.02);
  auto f=g; f.addIa(1.);
  auto old_f=old; old_f.addIa(1.);
  const auto acceleration=((g-old)/(dt*dt)-v/dt-a*(.5-beta))/beta;
  const auto velocity=v+a*(dt*(1-gamma))+acceleration*(dt*gamma);
  const auto l=velocity*f.inverse(), lo=v*old_f.inverse();
  const auto d=(l+l.transpose())*.5, d_old=(lo+lo.transpose())*.5;
  const Real mu=E/(2*(1+nu)), lambda=E*nu/((1+nu)*(1-2*nu));
  auto s=d*(2*mu); s.addIa(lambda*d.trace());
  auto so=d_old*(2*mu); so.addIa(lambda*d_old.trace());
  const auto expected=(((s*(1+alpha)-so*alpha)*zeta)*f.inverse().transpose())*f.det();
  const auto actual=FinitePlaneStrainVelocityRate::nominalStress(
      g,old,v,a,dt,beta,gamma,alpha,zeta,E,nu);
  EXPECT_LE((actual-expected).L2norm(),1e-12*expected.L2norm());
}

TEST(FinitePlaneStrainVelocityRate, FullFinitePullbackADMatchesFourGradientPerturbations)
{
  const auto g=tensor(.12,.03,-.02,-.04), old=tensor(.119,.029,-.019,-.039);
  const auto v=tensor(.004,.002,-.003,.001), a=tensor(.01,.02,.03,-.02);
  ADRankTwoTensor ad_g=g;
  for (unsigned int k=0;k<4;++k)
    Moose::derivInsert(ad_g(k/2,k%2).derivatives(),k,1.);
  const auto ad_p=FinitePlaneStrainVelocityRate::nominalStress(
      ad_g,old,v,a,dt,beta,gamma,alpha,zeta,E,nu);
  for (unsigned int k=0;k<4;++k)
  {
    auto plus=g, minus=g;
    constexpr Real h=1e-7;
    plus(k/2,k%2)+=h; minus(k/2,k%2)-=h;
    const auto p=FinitePlaneStrainVelocityRate::nominalStress(
        plus,old,v,a,dt,beta,gamma,alpha,zeta,E,nu);
    const auto m=FinitePlaneStrainVelocityRate::nominalStress(
        minus,old,v,a,dt,beta,gamma,alpha,zeta,E,nu);
    for (unsigned int i=0;i<2;++i)
      for (unsigned int j=0;j<2;++j)
      {
        const Real fd=(p(i,j)-m(i,j))/(2*h);
        EXPECT_NEAR(ad_p(i,j).derivatives()[k],fd,2e-7*std::max(1.,std::abs(fd)));
      }
  }
}

TEST(FinitePlaneStrainVelocityRate, ZeroDampingIsExactlyZeroIncludingDerivative)
{
  ADRankTwoTensor g=tensor(.12,.03,-.02,-.04);
  Moose::derivInsert(g(0,0).derivatives(),0,1.);
  const auto p=FinitePlaneStrainVelocityRate::nominalStress(
      g,RankTwoTensor(),RankTwoTensor(),RankTwoTensor(),dt,beta,gamma,alpha,0.,E,nu);
  EXPECT_DOUBLE_EQ(MetaPhysicL::raw_value(p.L2norm()),0.);
  EXPECT_DOUBLE_EQ(p(0,0).derivatives()[0],0.);
}

TEST(FinitePlaneStrainVelocityRate, RejectsInvertedCurrentAndOldConfigurations)
{
  const auto bad=tensor(-2.,0.,0.,0.), good=RankTwoTensor();
  EXPECT_THROW(FinitePlaneStrainVelocityRate::nominalStress(
      bad,good,good,good,dt,beta,gamma,alpha,zeta,E,nu),std::invalid_argument);
  EXPECT_THROW(FinitePlaneStrainVelocityRate::nominalStress(
      good,bad,good,good,dt,beta,gamma,alpha,zeta,E,nu),std::invalid_argument);
}
