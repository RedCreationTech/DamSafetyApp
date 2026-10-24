#include "gtest/gtest.h"
#include "AqueductElasticSlip.h"
#include "ADReal.h"

TEST(AqueductElasticSlip, ElasticRecoveryAndChangingPressure)
{
  auto a = AqueductFriction::update<double>({0, 0}, {.0003, .0004}, 100., .55, .001);
  EXPECT_NEAR(a.traction[0], 16.5, 1e-12);
  EXPECT_NEAR(a.traction[1], 22., 1e-12);
  auto b = AqueductFriction::update<double>(a.elastic, {0, 0}, 200., .55, .001);
  EXPECT_NEAR(b.traction[0], 33., 1e-12);
  auto c = AqueductFriction::update<double>(b.elastic, {-.0003, -.0004}, 200., .55, .001);
  EXPECT_NEAR(c.traction[0], 0., 1e-12);
}
TEST(AqueductElasticSlip, CircularReturnReverseAndOpenReset)
{
  auto a = AqueductFriction::update<double>({0, 0}, {.003, .004}, 100., .55, .001);
  EXPECT_NEAR(a.elastic[0], .0006, 1e-15);
  EXPECT_NEAR(a.elastic[1], .0008, 1e-15);
  EXPECT_NEAR(std::hypot(a.traction[0], a.traction[1]), 55., 1e-12);
  auto b = AqueductFriction::update<double>(a.elastic, {-.0006, -.0008}, 100., .55, .001);
  EXPECT_NEAR(std::hypot(b.traction[0], b.traction[1]), 0., 1e-12);
  auto c = AqueductFriction::update<double>(a.elastic, {9., 8.}, -1., .55, .001);
  EXPECT_EQ(c.elastic[0], 0.);
  EXPECT_EQ(c.traction[1], 0.);
}
TEST(AqueductElasticSlip, ConsistentPressureAndTwoDirectionDerivatives)
{
  for (double x : {.0002, .002})
  {
    ADReal a=x,b=.0003,p=100.;
    Moose::derivInsert(a.derivatives(), 0, 1.);
    Moose::derivInsert(b.derivatives(), 1, 1.);
    Moose::derivInsert(p.derivatives(), 2, 1.);
    auto state=AqueductFriction::update<ADReal>({0,0},{a,b},p,.55,.001);
    for (unsigned j=0;j<3;++j)
    {
      std::array<double,3> plus{x,.0003,100.},minus=plus;
      const double h=j==2?1e-3:1e-9;plus[j]+=h;minus[j]-=h;
      auto sp=AqueductFriction::update<double>({0,0},{plus[0],plus[1]},plus[2],.55,.001);
      auto sm=AqueductFriction::update<double>({0,0},{minus[0],minus[1]},minus[2],.55,.001);
      for (unsigned k=0;k<2;++k)
        EXPECT_NEAR(state.traction[k].derivatives()[j],(sp.traction[k]-sm.traction[k])/(2*h),1e-5);
    }
  }
}
