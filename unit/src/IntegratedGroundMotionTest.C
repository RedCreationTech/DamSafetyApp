#include "IntegratedGroundMotion.h"
#include "gtest/gtest.h"
#include <limits>

TEST(IntegratedGroundMotion, ConstantAccelerationWithInitialState)
{
  const IntegratedGroundMotion motion({0,1,3},{2,2,2},4,5);
  for (double t : {0.0,0.37,1.0,1.8,3.0})
  {
    const auto s=motion.at(t);
    EXPECT_NEAR(s.displacement,4+5*t+t*t,1e-13);
    EXPECT_NEAR(s.velocity,5+2*t,1e-13);
    EXPECT_DOUBLE_EQ(s.acceleration,2);
  }
}
TEST(IntegratedGroundMotion, LinearAccelerationIndependentOfSamplingPartition)
{
  const IntegratedGroundMotion coarse({0,2},{-0.004,1.996});
  const IntegratedGroundMotion fine({0,.07,.13,.51,1,2},{-.004,.066,.126,.506,.996,1.996});
  for(double t : {0.0,.013,.07,.131,.333,1.0,1.719,2.0})
  {
    const auto a=coarse.at(t), b=fine.at(t);
    EXPECT_NEAR(a.displacement,-.002*t*t+t*t*t/6,1e-14);
    EXPECT_NEAR(a.velocity,-.004*t+t*t/2,1e-14);
    EXPECT_NEAR(a.displacement,b.displacement,1e-14);
    EXPECT_NEAR(a.velocity,b.velocity,1e-14);
    EXPECT_NEAR(a.acceleration,b.acceleration,1e-14);
  }
}
TEST(IntegratedGroundMotion, CrossKnotTriangleHasKnownIntegral)
{
  const IntegratedGroundMotion motion({0,1,2},{0,1,0});
  const auto s=motion.at(2);
  EXPECT_NEAR(s.velocity,1,1e-14);
  EXPECT_NEAR(s.displacement,1,1e-14);
  EXPECT_NEAR(motion.at(1).displacement,1.0/6,1e-14);
  EXPECT_NEAR(motion.at(1).velocity,.5,1e-14);
}
TEST(IntegratedGroundMotion, QueryOrderAndRejectedTrialCannotChangeMotion)
{
  const IntegratedGroundMotion motion({0,.07,.13,.2},{-.004,.2,-.1,.05});
  const auto before=motion.at(.08);
  for(double t : {.2,.1,.09,.13,0.0,.2,.01}) (void)motion.at(t);
  const auto after=motion.at(.08);
  EXPECT_DOUBLE_EQ(before.displacement,after.displacement);
  EXPECT_DOUBLE_EQ(before.velocity,after.velocity);
  EXPECT_DOUBLE_EQ(before.acceleration,after.acceleration);
}
TEST(IntegratedGroundMotion, RejectsMalformedInputAndOutsideDomain)
{
  EXPECT_THROW(IntegratedGroundMotion({0},{1}),std::invalid_argument);
  EXPECT_THROW(IntegratedGroundMotion({0,0},{1,2}),std::invalid_argument);
  EXPECT_THROW(IntegratedGroundMotion({1,0},{1,2}),std::invalid_argument);
  EXPECT_THROW(IntegratedGroundMotion({0,1},{1}),std::invalid_argument);
  EXPECT_THROW(IntegratedGroundMotion({0,1},{1,std::numeric_limits<double>::infinity()}),std::invalid_argument);
  const IntegratedGroundMotion motion({0,1},{1,2});
  EXPECT_THROW(motion.at(-.001),std::out_of_range);
  EXPECT_THROW(motion.at(1.001),std::out_of_range);
}
