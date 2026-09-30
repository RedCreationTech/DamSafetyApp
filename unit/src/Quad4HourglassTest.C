#include "Quad4Hourglass.h"

#include "gtest/gtest.h"

namespace
{
double
dot(const Quad4Hourglass::Mode & left, const std::array<double, 4> & right)
{
  double value = 0.0;
  for (unsigned int i = 0; i < 4; ++i)
    value += left[i] * right[i];
  return value;
}
}

TEST(Quad4Hourglass, SquareMode)
{
  const Quad4Hourglass::Nodes nodes = {{{{0, 0}}, {{2, 0}}, {{2, 2}}, {{0, 2}}}};
  const auto gamma = Quad4Hourglass::mode(nodes);
  EXPECT_NEAR(std::abs(gamma[0]), 0.5, 1e-14);
  EXPECT_NEAR(dot(gamma, {{1, 1, 1, 1}}), 0.0, 1e-14);
  EXPECT_NEAR(dot(gamma, {{0, 2, 2, 0}}), 0.0, 1e-14);
  EXPECT_NEAR(dot(gamma, {{0, 0, 2, 2}}), 0.0, 1e-14);
  EXPECT_GT(0.5 * Quad4Hourglass::area(nodes) /
                Quad4Hourglass::characteristicLengthSquared(nodes) * dot(gamma, gamma),
            0.0);
  EXPECT_NEAR(Quad4Hourglass::area(nodes), 4.0, 1e-14);
  EXPECT_NEAR(Quad4Hourglass::characteristicLengthSquared(nodes), 4.0, 1e-14);
}

TEST(Quad4Hourglass, DistortedElementRejectsAllAffineFields)
{
  const Quad4Hourglass::Nodes nodes = {{{{0.1, -0.2}}, {{2.2, 0.1}}, {{1.8, 1.7}}, {{-0.3, 1.2}}}};
  const auto gamma = Quad4Hourglass::mode(nodes);
  std::array<double, 4> constant;
  std::array<double, 4> affine;
  for (unsigned int i = 0; i < 4; ++i)
  {
    constant[i] = 1.0;
    affine[i] = 3.0 - 2.0 * nodes[i][0] + 0.7 * nodes[i][1];
  }
  EXPECT_NEAR(dot(gamma, constant), 0.0, 1e-13);
  EXPECT_NEAR(dot(gamma, affine), 0.0, 1e-13);
  EXPECT_NEAR(dot(gamma, gamma), 1.0, 1e-13);
}

TEST(Quad4Hourglass, RejectsDegenerateElement)
{
  const Quad4Hourglass::Nodes nodes = {{{{0, 0}}, {{1, 0}}, {{2, 0}}, {{3, 0}}}};
  EXPECT_THROW(Quad4Hourglass::mode(nodes), std::runtime_error);
  EXPECT_THROW(Quad4Hourglass::area(nodes), std::runtime_error);
}

TEST(Quad4Hourglass, RestoresExactConsistentMassForSquare)
{
  const Quad4Hourglass::Nodes nodes = {{{{0, 0}}, {{2, 0}}, {{2, 2}}, {{0, 2}}}};
  const auto correction = Quad4Hourglass::consistentMassCorrection(nodes);
  for (unsigned int i = 0; i < 4; ++i)
  {
    double row_sum = 0.0;
    for (const auto value : correction[i])
      row_sum += value;
    EXPECT_NEAR(row_sum, 0.0, 1e-14);
    EXPECT_NEAR(correction[i][i], 7.0 / 36.0, 1e-14);
    EXPECT_NEAR(correction[i][(i + 1) % 4], -1.0 / 36.0, 1e-14);
    EXPECT_NEAR(correction[i][(i + 2) % 4], -5.0 / 36.0, 1e-14);
  }
}
