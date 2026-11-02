#include "CDPObservationFingerprint.h"
#include "gtest/gtest.h"

TEST(DS04GObservationFingerprint, TracksOwnedRangeOrderAndBits)
{
  auto run = [](std::uint64_t first, std::uint64_t last, double a, double b)
  {
    CDPObservationFingerprint h(first, last);
    h.add(a);
    h.add(b);
    return h.value();
  };
  const auto baseline = run(20, 22, 1.0, 0.0);
  EXPECT_EQ(baseline, run(20, 22, 1.0, 0.0));
  EXPECT_NE(baseline, run(20, 22, 0.0, 1.0));
  EXPECT_NE(baseline, run(21, 23, 1.0, 0.0));
  EXPECT_NE(baseline, run(20, 22, 1.0, -0.0));
}
