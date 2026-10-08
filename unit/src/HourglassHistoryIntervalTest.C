#include "StressHistoryInterval.h"
#include "Quad4Hourglass.h"
#include "gtest/gtest.h"

TEST(HourglassHistoryInterval, ConstantModalRateUsesBothActualIntervals)
{
  const double alpha = -.05, zeta = .03, rate = .7;
  for (double dt : {.01, .005, .0001})
    for (double previous_dt : {.02, .004, .00005})
    {
      const double current = .3, old = current - rate * dt, older = old - rate * previous_dt;
      const double original = (1 + alpha) * current - alpha * old +
          zeta * ((1 + alpha) * (current - old) / dt - alpha * (old - older) / dt);
      const double paired = original + StressHistoryInterval::correction(alpha, zeta, dt, previous_dt, 2) *
          (older - old);
      EXPECT_NEAR(paired, (1 + alpha) * current - alpha * old + zeta * rate, 2e-14);
    }
}

TEST(HourglassHistoryInterval, DoesNotChangeCurrentModalTangent)
{
  const double dt = .005, previous_dt = .02, alpha = -.05, zeta = .03;
  auto residual = [&](double current) {
    return (1 + alpha) * current - alpha * .2 + zeta *
        ((1 + alpha) * (current - .2) / dt - alpha * (.2 - .1) / dt) +
        StressHistoryInterval::correction(alpha, zeta, dt, previous_dt, 3) * (.1 - .2);
  };
  EXPECT_NEAR((residual(.301) - residual(.299)) / .002, (1 + alpha) * (1 + zeta / dt), 1e-12);
}

TEST(HourglassHistoryInterval, AffineHistoryRemainsUnpenalized)
{
  const Quad4Hourglass::Nodes nodes = {{{{0, 0}}, {{2, 0}}, {{1, 1}}, {{0, 1}}}};
  const auto gamma = Quad4Hourglass::mode(nodes);
  double old_mode = 0, older_mode = 0;
  for (unsigned int i = 0; i < 4; ++i)
  {
    old_mode += gamma[i] * (3 + .7 * nodes[i][0] - .9 * nodes[i][1]);
    older_mode += gamma[i] * (2 + .4 * nodes[i][0] - .6 * nodes[i][1]);
  }
  EXPECT_NEAR(StressHistoryInterval::correction(-.05, .03, .005, .01, 3) *
              (older_mode - old_mode), 0, 1e-14);
}
