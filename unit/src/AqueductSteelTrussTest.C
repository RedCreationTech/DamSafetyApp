#include "AqueductSteelTruss.h"
#include "gtest/gtest.h"

namespace {constexpr double E=200e9,fy=400e6,fu=696e6,limit=.1;}
TEST(AqueductSteelTruss, MonotonicClosedFormAndConsistentTangent)
{
  auto s=AqueductTruss::update(.01,0,0,E,fy,fu,limit);
  double H=(fu-fy)/limit,p=(E*.01-fy)/(E+H);
  EXPECT_NEAR(s.stress,fy+H*p,1e-5);EXPECT_NEAR(s.equivalent,p,1e-14);
  auto a=AqueductTruss::update(.01+1e-7,0,0,E,fy,fu,limit);
  auto b=AqueductTruss::update(.01-1e-7,0,0,E,fy,fu,limit);
  EXPECT_NEAR(s.tangent,(a.stress-b.stress)/2e-7,1);
}
TEST(AqueductSteelTruss, UnloadAndReverseAccumulateIsotropicHardening)
{
  auto s=AqueductTruss::update(.01,0,0,E,fy,fu,limit);
  auto u=AqueductTruss::update(.0099,s.plastic,s.equivalent,E,fy,fu,limit);
  EXPECT_NEAR(u.stress,s.stress-E*.0001,1e-5);EXPECT_DOUBLE_EQ(u.equivalent,s.equivalent);
  auto r=AqueductTruss::update(-.01,u.plastic,u.equivalent,E,fy,fu,limit);
  EXPECT_GT(r.equivalent,u.equivalent);
  EXPECT_NEAR(std::abs(r.stress),fy+(fu-fy)/limit*r.equivalent,1e-5);
  EXPECT_GT(r.equivalent,std::abs(r.plastic));
}
TEST(AqueductSteelTruss, SourceLastRowUsesConstantExtrapolation)
{
  auto s=AqueductTruss::update(.2,0,0,E,fy,fu,limit);
  EXPECT_NEAR(s.stress,fu,1e-5);EXPECT_DOUBLE_EQ(s.tangent,0);
  auto r=AqueductTruss::update(-.2,s.plastic,s.equivalent,E,fy,fu,limit);
  EXPECT_NEAR(r.stress,-fu,1e-5);EXPECT_GT(r.equivalent,s.equivalent);
}
