#include "UniformHex8.h"
#include "gtest/gtest.h"
#include <stdexcept>

namespace
{
std::array<libMesh::Point, 8> cube()
{
  return {{{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}}};
}
}

TEST(UniformHex8, CubeVolumeAndConstantGradient)
{
  const auto op = UniformHex8::compute(cube());
  EXPECT_NEAR(op.volume, 1.0, 1e-14);
  for (unsigned int d = 0; d < 3; ++d)
  {
    Real sum = 0;
    for (const auto & b : op.gradient) sum += b(d);
    EXPECT_NEAR(sum, 0.0, 1e-14);
  }
}

TEST(UniformHex8, DistortedBrickPassesAffinePatchAndHourglassNullTest)
{
  auto points = cube();
  points[6] += libMesh::Point(0.17,-0.08,0.12);
  for (auto & p : points) p += libMesh::Point(3.2,-8.1,6.3);
  const auto op = UniformHex8::compute(points);
  for (unsigned int a = 0; a < 3; ++a)
    for (unsigned int b = 0; b < 3; ++b)
    {
      Real gradient = 0;
      for (unsigned int i = 0; i < 8; ++i) gradient += points[i](a) * op.gradient[i](b);
      EXPECT_NEAR(gradient, a == b ? 1.0 : 0.0, 3e-14);
    }
  for (unsigned int i = 0; i < 8; ++i)
    for (unsigned int d = 0; d < 4; ++d)
    {
      Real force = 0;
      for (unsigned int j = 0; j < 8; ++j)
        force += op.hourglassStiffness(i,j,4e8) * (d == 3 ? 1.0 : points[j](d));
      EXPECT_NEAR(force, 0.0, 2e-7);
    }
}

TEST(UniformHex8, FourHourglassModesHavePositiveEnergy)
{
  const auto op = UniformHex8::compute(cube());
  for (unsigned int mode = 0; mode < 4; ++mode)
  {
    Real energy = 0;
    for (unsigned int i = 0; i < 8; ++i)
      for (unsigned int j = 0; j < 8; ++j)
      {
        EXPECT_NEAR(op.hourglassStiffness(i,j,4e8),op.hourglassStiffness(j,i,4e8),1e-8);
        energy += op.hourglass[i][mode] * op.hourglassStiffness(i,j,4e8) * op.hourglass[j][mode];
      }
    EXPECT_GT(energy, 0.0);
  }
}

TEST(UniformHex8, InvalidOrientationFails)
{
  auto points = cube();
  for (auto & p : points) p(0) = -p(0);
  EXPECT_THROW(UniformHex8::compute(points), std::invalid_argument);
}
