#include "AqueductOutputSchedule.h"
#include "gtest/gtest.h"
#include <limits>

TEST(AqueductOutputSchedule, MarkAnchoringDoesNotDriftAfterCutbacks)
{
  auto slot=[](double t){return AqueductOutputSchedule::slot(t,2.,12.,.1,1e-10);};
  EXPECT_EQ(slot(1.99),-1);
  EXPECT_EQ(slot(2.),0);
  EXPECT_EQ(slot(2.099),0);
  EXPECT_EQ(slot(2.1013),1);
  EXPECT_EQ(slot(2.2001),2); // Unlike last-written-time + interval, no accumulating delay.
  EXPECT_EQ(slot(12.),100);
  EXPECT_EQ(slot(12.001),-1);
}
TEST(AqueductOutputSchedule, InvalidScheduleIsRejected)
{
  EXPECT_THROW(AqueductOutputSchedule::slot(2.,2.,12.,0.,1e-10),std::invalid_argument);
  EXPECT_THROW(AqueductOutputSchedule::slot(2.,2.,1.,.1,1e-10),std::invalid_argument);
  EXPECT_THROW(AqueductOutputSchedule::slot(2.,2.,12.,.1,.1),std::invalid_argument);
  EXPECT_THROW(AqueductOutputSchedule::slot(std::numeric_limits<double>::quiet_NaN(),2.,12.,.1,1e-10),std::invalid_argument);
}
