#include "StressHistoryInterval.h"
#include "gtest/gtest.h"

TEST(StressHistoryInterval, EqualIntervalsRetainUpstreamFormula)
{
  EXPECT_DOUBLE_EQ(StressHistoryInterval::correction(-.05,.00113,.01,.01,2),0);
  EXPECT_DOUBLE_EQ(StressHistoryInterval::correction(-.05,.00113,.01,0,1),0);
  EXPECT_DOUBLE_EQ(StressHistoryInterval::correction(-.05,.00113,0,0,0),0);
}
TEST(StressHistoryInterval, ConstantRateIsIndependentOfPreviousInterval)
{
  const double alpha=-.05,zeta=.00113,rate=713;
  for (double dt : {.01,.005,.0001})
    for (double old_dt : {.02,.005,.00005})
    {
      const double stress=17,old=stress-rate*dt,older=old-rate*old_dt;
      const double upstream=(1+alpha)*stress-alpha*old+
          zeta*((1+alpha)*(stress-old)/dt-alpha*(old-older)/dt);
      const double fixed=upstream+StressHistoryInterval::correction(alpha,zeta,dt,old_dt,2)*(older-old);
      EXPECT_NEAR(fixed,(1+alpha)*stress-alpha*old+zeta*rate,1e-11);
    }
}
TEST(StressHistoryInterval, HistoryCorrectionDoesNotChangeCurrentTangent)
{
  const double c=StressHistoryInterval::correction(-.05,.03,.01,.025,3);
  const double old=4,older=1;
  auto residual=[&](double current) { return 3.8*current+c*(older-old); };
  EXPECT_NEAR((residual(2.001)-residual(1.999))/.002,3.8,1e-12);
}
TEST(StressHistoryInterval, InvalidAcceptedHistoryIsRejected)
{
  EXPECT_THROW(StressHistoryInterval::correction(-.05,.03,.01,0,2),std::invalid_argument);
  EXPECT_DOUBLE_EQ(StressHistoryInterval::correction(0,.03,.01,0,2),0);
}
