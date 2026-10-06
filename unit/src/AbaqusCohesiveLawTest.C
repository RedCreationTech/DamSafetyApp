#include "AbaqusCohesiveLaw.h"
#include "gtest/gtest.h"

namespace
{
const AbaqusCohesiveLaw::Parameters law{{2.559e12,1.093e12,1.093e12},
                                      {5.63e6,5.63e6,5.63e6},0.000241,0,0.6,0.0003};
}
TEST(AbaqusCohesiveLaw, IndependentNormalAndShearStiffnessAndCompression)
{
  auto out=AbaqusCohesiveLaw::evaluate<double>(law,{},{{1e-7,2e-7,-3e-7}},0.0,0.01);
  EXPECT_DOUBLE_EQ(out.damage,0);
  EXPECT_NEAR(out.traction[0],-255900,1e-8);
  EXPECT_NEAR(out.traction[1],-218600,1e-8);
  EXPECT_NEAR(out.traction[2],327900,1e-8);
  out=AbaqusCohesiveLaw::evaluate<double>(law,{},{{-0.01,0,0}},100.0,0.01);
  EXPECT_DOUBLE_EQ(out.damage,0);EXPECT_DOUBLE_EQ(out.traction[0],100);
}
TEST(AbaqusCohesiveLaw, MaximumStressIsNotQuadraticCriterion)
{
  const double n=0.9*law.strength[0]/law.stiffness[0];
  const double s=0.9*law.strength[1]/law.stiffness[1];
  const auto out=AbaqusCohesiveLaw::evaluate<double>(law,{},{{n,s,s}},0.0,0.01);
  EXPECT_DOUBLE_EQ(out.state[0],0);EXPECT_DOUBLE_EQ(out.damage,0);
}
TEST(AbaqusCohesiveLaw, LinearSofteningUsesFailureIncrementAfterOnset)
{
  const double onset=law.strength[0]/law.stiffness[0];
  auto first=AbaqusCohesiveLaw::evaluate<double>(law,{},{{onset,0,0}},0.0,0.01);
  EXPECT_NEAR(first.state[0],onset,1e-15);EXPECT_NEAR(first.traction[0],-law.strength[0],1e-6);
  const auto middle=AbaqusCohesiveLaw::evaluate<double>(law,first.state,{{onset+0.5*law.failure_increment,0,0}},0.0,0.01);
  EXPECT_NEAR(middle.traction[0],-0.5*law.strength[0],1e-6);
  const auto end=AbaqusCohesiveLaw::evaluate<double>(law,middle.state,{{onset+law.failure_increment,0,0}},0.0,0.01);
  EXPECT_NEAR(end.damage,1,1e-14);EXPECT_NEAR(end.traction[0],0,1e-6);
  const auto unload=AbaqusCohesiveLaw::evaluate<double>(law,middle.state,{{0.25*onset,0,0}},0.0,0.01);
  EXPECT_NEAR(unload.damage,middle.damage,1e-14);
}
TEST(AbaqusCohesiveLaw, RejectedTrialDoesNotCommitInitiationOrDamage)
{
  const AbaqusCohesiveLaw::State old{};
  const auto large=AbaqusCohesiveLaw::evaluate<double>(law,old,{{0.01,0,0}},0.0,0.01);
  EXPECT_DOUBLE_EQ(large.damage,1);
  const auto retry=AbaqusCohesiveLaw::evaluate<double>(law,old,{{1e-7,0,0}},0.0,0.005);
  EXPECT_DOUBLE_EQ(retry.damage,0);EXPECT_DOUBLE_EQ(retry.state[0],0);
}
TEST(AbaqusCohesiveLaw, ViscosityAndDamageWeightedCoulombFriction)
{
  auto viscous=law;viscous.viscosity=0.0005;
  auto old=AbaqusCohesiveLaw::State{};
  old[0]=law.strength[0]/law.stiffness[0];old[1]=old[0]+law.failure_increment;old[8]=1;
  const auto out=AbaqusCohesiveLaw::evaluate<double>(viscous,old,{{-1e-5,0.001,0}},1000.0,0.0005);
  EXPECT_NEAR(out.damage,0.5,1e-14);
  EXPECT_NEAR(out.traction[1],-0.5*law.stiffness[1]*0.001-0.5*600,1e-7);
  const auto failed=AbaqusCohesiveLaw::evaluate<double>(law,{},{{0,0.001,0}},1000.0,0.01,false);
  EXPECT_NEAR(failed.traction[1],-600,1e-10);
}
TEST(AbaqusCohesiveLaw, AutomaticDerivativeMatchesFiniteDifferenceOnSofteningBranch)
{
  const double displacement=law.strength[0]/law.stiffness[0]+0.4*law.failure_increment;
  std::array<ADReal,3> jump{displacement,0,0};Moose::derivInsert(jump[0].derivatives(),0,1.0);
  const auto ad=AbaqusCohesiveLaw::evaluate<ADReal>(law,{},jump,ADReal(0),0.01);
  const double h=1e-9;
  const auto a=AbaqusCohesiveLaw::evaluate<double>(law,{},{{displacement+h,0,0}},0.0,0.01);
  const auto b=AbaqusCohesiveLaw::evaluate<double>(law,{},{{displacement-h,0,0}},0.0,0.01);
  EXPECT_NEAR(ad.traction[0].derivatives()[0],(a.traction[0]-b.traction[0])/(2*h),1e3);
}
