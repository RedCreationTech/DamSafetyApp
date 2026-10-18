#include "Quad4Hourglass.h"
#include "StressHistoryInterval.h"
#include "MooseTypes.h"
#include "ADUtils.h"
#include "gtest/gtest.h"
#include <algorithm>

namespace
{
const Quad4Hourglass::Nodes square = {{{{0, 0}}, {{1, 0}}, {{1, 1}}, {{0, 1}}}};
const Quad4Hourglass::Nodes distorted = {{{{0.1, -0.2}}, {{2.2, 0.1}},
                                         {{1.8, 1.7}}, {{-0.3, 1.2}}}};
ADReal seeded(double value, unsigned int index)
{
  ADReal x = value;
  Moose::derivInsert(x.derivatives(), index, 1.0);
  return x;
}
}

TEST(HourglassEnergyReference, SquareAndRectangleAbsoluteCoefficient)
{
  EXPECT_NEAR(Quad4Hourglass::gradientRawStiffnessPerShear(square), 0.0025, 1e-17);
  const Quad4Hourglass::Nodes rectangle = {{{{0, 0}}, {{2, 0}}, {{2, 1}}, {{0, 1}}}};
  EXPECT_NEAR(Quad4Hourglass::gradientRawStiffnessPerShear(rectangle), 0.0015625, 1e-17);
}

TEST(HourglassEnergyReference, ScaleRotationAndTranslation)
{
  auto transformed = distorted;
  const double angle = 0.7, scale = 3.7;
  for (unsigned int i = 0; i < 4; ++i)
  {
    transformed[i][0] = scale * (std::cos(angle)*distorted[i][0] -
                                  std::sin(angle)*distorted[i][1]) + 21.;
    transformed[i][1] = scale * (std::sin(angle)*distorted[i][0] +
                                  std::cos(angle)*distorted[i][1]) + 75.;
  }
  EXPECT_NEAR(Quad4Hourglass::gradientRawStiffnessPerShear(transformed)*scale*scale,
              Quad4Hourglass::gradientRawStiffnessPerShear(distorted), 1e-16);
}

TEST(HourglassEnergyReference, ActualADNodalMatrixAndDirectionalSlopes)
{
  const auto gamma = Quad4Hourglass::mode(distorted);
  ADReal q = 0.0;
  for (unsigned int j = 0; j < 4; ++j)
    q += gamma[j] * seeded(.001*(j+1), j);
  const double alpha = -.05, zeta = .00113, G = 17., old = -.002, older = .001;
  const double kref = Quad4Hourglass::gradientRawStiffnessPerShear(distorted);
  const double kold = 1./Quad4Hourglass::characteristicLengthSquared(distorted);
  for (double dt : {0., .01, .005, .000234375, .0000005})
    for (double previous_dt : {.01, .003, .02})
    {
      const double correction = dt > 0 ? StressHistoryInterval::correction(
          alpha, zeta, dt, previous_dt, 3) : 0.;
      const double slope = dt > 0 ? (1+alpha)*(kref+kold*zeta/dt) : kref;
      for (unsigned int i = 0; i < 4; ++i)
      {
        const ADReal force = G*gamma[i]*Quad4Hourglass::gradientRawResidualPerShear(
            distorted, alpha, dt, zeta, correction, q, old, older);
        for (unsigned int j = 0; j < 4; ++j)
          EXPECT_NEAR(force.derivatives()[j], G*gamma[i]*gamma[j]*slope,
                      1e-11*std::max(1.,std::abs(force.derivatives()[j])));
      }
      auto residual = [&](double value) {
        return Quad4Hourglass::gradientRawResidualPerShear(
            distorted, alpha, dt, zeta, correction, value, old, older);
      };
      const double value = .003, eps = 1e-6;
      EXPECT_NEAR((residual(value+eps)-residual(value))/eps, slope,
                  1e-9*std::max(1.,std::abs(slope)));
      EXPECT_NEAR((residual(value)-residual(value-eps))/eps, slope,
                  1e-9*std::max(1.,std::abs(slope)));
    }
}

TEST(HourglassEnergyReference, EntireLegacyDampingValueAndADAreRetained)
{
  ADReal q = seeded(.003, 0);
  const double alpha = -.05, zeta = .00113, old = -.002, older = .001;
  const double kold = 1./Quad4Hourglass::characteristicLengthSquared(distorted);
  for (double dt : {.01, .005, .000234375, .0000005})
    for (double previous_dt : {.01, .003, .02})
      for (bool use_previous : {false, true})
      {
        const double correction = use_previous ? StressHistoryInterval::correction(
            alpha, zeta, dt, previous_dt, 3) : 0.;
        const ADReal original = kold*(q*((1+alpha)*zeta/dt) -
            old*((1+2*alpha)*zeta/dt) + older*(alpha*zeta/dt) +
            correction*(older-old));
        const ADReal on = Quad4Hourglass::gradientRawResidualPerShear(
            distorted, alpha, dt, zeta, correction, q, old, older);
        const ADReal off = Quad4Hourglass::gradientRawResidualPerShear(
            distorted, alpha, dt, 0., 0., q, old, older);
        EXPECT_NEAR((on-off).value(), original.value(), 1e-12);
        EXPECT_NEAR((on-off).derivatives()[0], original.derivatives()[0], 1e-10);
      }
}

TEST(HourglassEnergyReference, StaticLimitAndAffineNullEnergy)
{
  ADReal q = seeded(.003, 0);
  const double k = Quad4Hourglass::gradientRawStiffnessPerShear(distorted);
  const ADReal force = Quad4Hourglass::gradientRawResidualPerShear(
      distorted, -.05, 0., .03, 7., q, 17., -19.);
  EXPECT_NEAR(force.value(), k*.003, 1e-18);
  EXPECT_NEAR(force.derivatives()[0], k, 1e-18);
  const auto gamma = Quad4Hourglass::mode(distorted);
  double rigid = 0., affine = 0.;
  for (unsigned int i = 0; i < 4; ++i)
  {
    rigid += gamma[i]*7.;
    affine += gamma[i]*(3. + .7*distorted[i][0] - .9*distorted[i][1]);
  }
  EXPECT_NEAR(rigid, 0., 1e-13);
  EXPECT_NEAR(affine, 0., 1e-13);
  EXPECT_GT(k, 0.);
}

TEST(HourglassEnergyReference, RejectedTrialCannotModifyAcceptedHistory)
{
  const double old = -.002, older = .001;
  const auto accepted = std::array<double, 2>{{old, older}};
  const auto retry = [&]() {
    return Quad4Hourglass::gradientRawResidualPerShear(
        distorted, -.05, .005, .00113, 0., seeded(.003,0), old, older);
  };
  const ADReal before = retry();
  const ADReal rejected = Quad4Hourglass::gradientRawResidualPerShear(
      distorted, -.05, .01, .00113, 0., seeded(1.5,0), old, older);
  EXPECT_NE(rejected.value(), before.value());
  const ADReal after = retry();
  EXPECT_EQ(before.value(), after.value());
  EXPECT_EQ(before.derivatives()[0], after.derivatives()[0]);
  EXPECT_EQ(accepted[0], old);
  EXPECT_EQ(accepted[1], older);
}

TEST(HourglassEnergyReference, RejectsNonpositiveCenterJacobian)
{
  const Quad4Hourglass::Nodes flat = {{{{0,0}},{{1,0}},{{2,0}},{{3,0}}}};
  auto reversed = square;
  std::swap(reversed[1], reversed[3]);
  EXPECT_THROW(Quad4Hourglass::gradientRawStiffnessPerShear(flat), std::runtime_error);
  EXPECT_THROW(Quad4Hourglass::gradientRawStiffnessPerShear(reversed), std::runtime_error);
}
