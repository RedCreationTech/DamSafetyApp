#include "HourglassElasticScale.h"
#include "Quad4Hourglass.h"
#include "StressHistoryInterval.h"
#include "MooseTypes.h"
#include "ADUtils.h"
#include "gtest/gtest.h"

TEST(HourglassElasticScale, DefaultAddsNeitherForceNorADDerivative)
{
  ADReal current = .003;
  Moose::derivInsert(current.derivatives(), 0, 1.0);
  for (double dt : {0., .0000005, .01})
  {
    const ADReal result = HourglassElasticScale::correction(1., -.05, dt, current, -.001);
    EXPECT_EQ(result.value(), 0.);
    EXPECT_EQ(result.derivatives()[0], 0.);
  }
}

TEST(HourglassElasticScale, ModalMatrixHasOnlyTheElasticRankOneChange)
{
  const Quad4Hourglass::Nodes nodes = {{{{0., 0.}}, {{2., 0.}}, {{1., 1.}}, {{0., 1.}}}};
  const auto gamma = Quad4Hourglass::mode(nodes);
  ADReal mode = 0.;
  for (unsigned int j = 0; j < 4; ++j)
  {
    ADReal value = .001 * (j + 1);
    Moose::derivInsert(value.derivatives(), j, 1.);
    mode += gamma[j] * value;
  }
  for (unsigned int i = 0; i < 4; ++i)
  {
    const ADReal force = 17. * gamma[i] *
        HourglassElasticScale::correction(.8, -.05, .003, mode, -.001);
    for (unsigned int j = 0; j < 4; ++j)
      EXPECT_NEAR(force.derivatives()[j], 17. * (-.2) * .95 * gamma[i] * gamma[j], 1e-13);
  }
}

TEST(HourglassElasticScale, AdditionalDampingForceAndADDerivativeAreUnchanged)
{
  const double alpha = -.05, old = -.002, older = .001;
  ADReal current = .003;
  Moose::derivInsert(current.derivatives(), 0, 1.);
  for (double dt : {.0000005, .000234375, .005, .01})
    for (double previous_dt : {.00005, .003, .02})
      for (double scale : {.8, 1., 1.2})
      {
        auto residual = [&](double zeta, double s) -> ADReal {
          return current * (1. + alpha + (1. + alpha) * zeta / dt) -
              old * (alpha + (1. + 2. * alpha) * zeta / dt) + older * alpha * zeta / dt +
              StressHistoryInterval::correction(alpha, zeta, dt, previous_dt, 3) * (older - old) +
              HourglassElasticScale::correction(s, alpha, dt, current, old);
        };
        const ADReal original_damping = residual(.00113, 1.) - residual(0., 1.);
        const ADReal scaled_damping = residual(.00113, scale) - residual(0., scale);
        EXPECT_NEAR(original_damping.value(), scaled_damping.value(), 1e-12);
        EXPECT_NEAR(original_damping.derivatives()[0], scaled_damping.derivatives()[0], 1e-10);
      }
}

TEST(HourglassElasticScale, StaticLimitAndAffineMotion)
{
  ADReal current = .003;
  Moose::derivInsert(current.derivatives(), 0, 1.);
  const ADReal delta = HourglassElasticScale::correction(.8, -.05, 0., current, 0.);
  EXPECT_NEAR(delta.value(), -.0006, 1e-16);
  EXPECT_NEAR(delta.derivatives()[0], -.2, 1e-15);
  const Quad4Hourglass::Nodes nodes = {{{{0., 0.}}, {{2., 0.}}, {{1., 1.}}, {{0., 1.}}}};
  const auto gamma = Quad4Hourglass::mode(nodes);
  double mode = 0., old = 0.;
  for (unsigned int i = 0; i < 4; ++i)
  {
    mode += gamma[i] * (3. + .7 * nodes[i][0] - .9 * nodes[i][1]);
    old += gamma[i] * (2. + .4 * nodes[i][0] - .6 * nodes[i][1]);
  }
  EXPECT_NEAR(HourglassElasticScale::correction(.8, -.05, .01, mode, old), 0., 1e-14);
}
