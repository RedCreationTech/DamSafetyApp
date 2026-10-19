#include "Quad4Hourglass.h"
#include "MooseTypes.h"
#include "ADUtils.h"
#include "gtest/gtest.h"
#include <algorithm>

namespace
{
ADReal seededRate(double x, unsigned int j)
{
  ADReal value = x;
  Moose::derivInsert(value.derivatives(), j, 1.0);
  return value;
}
const Quad4Hourglass::Nodes distortedRate = {{{{.1,-.2}},{{2.2,.1}},{{1.8,1.7}},{{-.3,1.2}}}};
}

TEST(HourglassNewmarkRate, ConstantVelocityAndStaticGuard)
{
  for (double h : {.01, .0021875, .000234375})
  {
    const double old=.003, v=.12;
    EXPECT_NEAR(Quad4Hourglass::newmarkWeightedRate(old+h*v,old,v,0.,h,.275625,.55,-.05),v,2e-14);
  }
  const ADReal zero = Quad4Hourglass::newmarkWeightedRate(seededRate(1.,0),3.,4.,5.,0.,.275625,.55,-.05);
  EXPECT_EQ(zero.value(),0.);
  EXPECT_EQ(zero.derivatives()[0],0.);
  EXPECT_THROW(Quad4Hourglass::newmarkWeightedRate(1.,0.,0.,0.,.01,0.,.55,-.05),std::invalid_argument);
}

TEST(HourglassNewmarkRate, IndependentNewmarkUpdate)
{
  for (double h : {.01, .0021875, .000234375})
    for (double beta : {.25,.275625})
      for (double gamma : {.5,.55})
      {
        const double old=.003,vold=-.12,aold=.7,anew=-1.3,alpha=-.05;
        const double current=old+h*vold+h*h*((.5-beta)*aold+beta*anew);
        const double velocity=vold+h*((1-gamma)*aold+gamma*anew);
        EXPECT_NEAR(Quad4Hourglass::newmarkWeightedRate(current,old,vold,aold,h,beta,gamma,alpha),
                    (1+alpha)*velocity-alpha*vold,1e-13);
      }
}

TEST(HourglassNewmarkRate, ActualADEightDofMatrixAndBothSlopes)
{
  const auto mode=Quad4Hourglass::mode(distortedRate);
  const double G=17.,zeta=.00113,alpha=-.05,beta=.275625,gamma=.55;
  const double ke=Quad4Hourglass::gradientRawStiffnessPerShear(distortedRate);
  const double kd=1./Quad4Hourglass::characteristicLengthSquared(distortedRate);
  for (double h : {.01,.005,.0021875,.000234375})
    for (unsigned int axis=0;axis<2;++axis)
    {
      ADReal q=0.,rate=0.;double qold=0.;
      for (unsigned int j=0;j<4;++j)
      {
        ADReal u=seededRate(.001*(j+1),4*axis+j);
        q+=mode[j]*u;qold+=mode[j]*(-.002*(j+1));
        rate+=mode[j]*Quad4Hourglass::newmarkWeightedRate(u,-.002*(j+1),.01*j,-.02*j,h,beta,gamma,alpha);
      }
      for (unsigned int i=0;i<4;++i)
      {
        ADReal force=G*mode[i]*(ke*((1+alpha)*q-alpha*qold)+kd*zeta*rate);
        for (unsigned int j=0;j<8;++j)
        {
          const double expected=j/4==axis?G*mode[i]*mode[j%4]*(1+alpha)*(ke+kd*zeta*gamma/(beta*h)):0.;
          EXPECT_NEAR(force.derivatives()[j],expected,1e-11*std::max(1.,std::abs(expected)));
        }
      }
      const double u=.003,old=-.002,eps=1e-7;
      const auto f=[&](double v){return Quad4Hourglass::newmarkWeightedRate(v,old,.01,-.02,h,beta,gamma,alpha);};
      const double slope=(1+alpha)*gamma/(beta*h);
      EXPECT_NEAR((f(u+eps)-f(u))/eps,slope,1e-8*std::abs(slope));
      EXPECT_NEAR((f(u)-f(u-eps))/eps,slope,1e-8*std::abs(slope));
    }
}

TEST(HourglassNewmarkRate, RigidAndAffineNullRate)
{
  const auto mode=Quad4Hourglass::mode(distortedRate);
  double sum=0.;
  for (unsigned int j=0;j<4;++j)
  {
    const double x=distortedRate[j][0],y=distortedRate[j][1];
    sum+=mode[j]*Quad4Hourglass::newmarkWeightedRate(1.+.003*x-.002*y, .8+.002*x-.001*y,
                                                   .01+.003*x+.002*y,-.02+.001*x-.003*y,
                                                   .01,.275625,.55,-.05);
  }
  EXPECT_NEAR(sum,0.,1e-9);
}

TEST(HourglassNewmarkRate, RepeatedRejectedEvaluationCannotAdvanceHistory)
{
  const double old=.003,vold=.12,aold=-.5;
  const ADReal before=Quad4Hourglass::newmarkWeightedRate(seededRate(.004,0),old,vold,aold,.005,.275625,.55,-.05);
  for (double h : {.01,.005,.0025})
  {
    const auto failed=Quad4Hourglass::newmarkWeightedRate(seededRate(.5,0),old,vold,aold,h,.275625,.55,-.05);
    EXPECT_NE(failed.value(),before.value());
  }
  const ADReal after=Quad4Hourglass::newmarkWeightedRate(seededRate(.004,0),old,vold,aold,.005,.275625,.55,-.05);
  EXPECT_EQ(before.value(),after.value());
  EXPECT_EQ(before.derivatives()[0],after.derivatives()[0]);
  EXPECT_EQ(old,.003);EXPECT_EQ(vold,.12);EXPECT_EQ(aold,-.5);
}
