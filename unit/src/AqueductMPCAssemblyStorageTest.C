#include "AqueductMPCAssemblyStorage.h"
#include "gtest/gtest.h"
#include <array>
#include <cmath>

TEST(AqueductMPCAssemblyStorage, CachedJacobianMatchesIndependentResidualDerivative)
{
  const std::vector<double> c{1., -.3, -.7, .1};
  std::array<double, 4> scales{2., .2, 1., 3.};
  const double multiplier_scale = .7, compliance = 1./605000.;
  AqueductMPCAssemblyStorage storage(c);
  const auto & matrix = storage.jacobian(scales, multiplier_scale, compliance);
  const std::array<double, 5> x{.01, -.02, .03, .004, 4.};
  auto residual = [&](const std::array<double, 5> & q) {
    std::array<double, 5> r{};
    double gap = -compliance * q.back();
    for (std::size_t i = 0; i < c.size(); ++i)
    {
      r[i] = scales[i] * c[i] * q.back();
      gap += c[i] * q[i];
    }
    r.back() = multiplier_scale * gap;
    return r;
  };
  const double h = 1e-6;
  for (std::size_t column = 0; column < x.size(); ++column)
  {
    auto plus = x, minus = x;
    plus[column] += h;
    minus[column] -= h;
    auto rp = residual(plus), rm = residual(minus);
    for (std::size_t row = 0; row < x.size(); ++row)
      EXPECT_NEAR(matrix(row, column), (rp[row]-rm[row])/(2*h), 1e-9);
  }
  for (unsigned i = 0; i < 100; ++i)
    EXPECT_EQ(&storage.jacobian(scales, multiplier_scale, compliance), &matrix);
  EXPECT_EQ(storage.rebuilds(), 1u);
}

TEST(AqueductMPCAssemblyStorage, ScalingAndComplianceChangesCannotUseStaleMatrix)
{
  AqueductMPCAssemblyStorage storage({1., -.2, -.8});
  std::array<double, 3> scales{1., 2., 3.};
  storage.jacobian(scales, 4., 0.);
  scales[1] = 7.;
  EXPECT_DOUBLE_EQ(storage.jacobian(scales, 4., 0.)(1,3), -1.4);
  EXPECT_EQ(storage.rebuilds(), 2u);
  EXPECT_DOUBLE_EQ(storage.jacobian(scales, 9., 0.)(3,2), -.8*9.);
  EXPECT_EQ(storage.rebuilds(), 3u);
  EXPECT_DOUBLE_EQ(storage.jacobian(scales, 9., .25)(3,3), -2.25);
  EXPECT_EQ(storage.rebuilds(), 4u);
  // Retry a previous scaling: the cached matrix must be rebuilt, not reused by step number.
  scales = {1., 2., 3.};
  EXPECT_DOUBLE_EQ(storage.jacobian(scales, 4., 0.)(1,3), -.4);
  EXPECT_EQ(storage.rebuilds(), 5u);
}

TEST(AqueductMPCAssemblyStorage, ScratchDoesNotRetainOldResidualOrDofValues)
{
  AqueductMPCAssemblyStorage storage({1., -.5, -.5});
  auto & ids = storage.dofs();
  ids = {10, 11, 12, 13};
  const auto capacity = ids.capacity();
  auto & next = storage.dofs();
  EXPECT_TRUE(next.empty());
  next.insert(next.end(), {21, 20, 23, 22});
  EXPECT_EQ(next[0], 21u);
  EXPECT_EQ(next.capacity(), capacity);
  storage.residual()(2) = 123.;
  const auto & r = storage.residual();
  for (unsigned j = 0; j < 4; ++j)
    EXPECT_DOUBLE_EQ(r(j), 0.);
}

TEST(AqueductMPCAssemblyStorage, EmbeddedAndPhysicalCompliancePreserveZeroConvention)
{
  std::array<double, 2> scales{1., 1.};
  AqueductMPCAssemblyStorage embedded({1., -1.}, false);
  AqueductMPCAssemblyStorage physical({1., -1.});
  EXPECT_FALSE(std::signbit(embedded.jacobian(scales, 1., 0.)(2,2)));
  EXPECT_TRUE(std::signbit(physical.jacobian(scales, 1., 0.)(2,2)));
  EXPECT_THROW(physical.jacobian(std::array<double, 1>{1.}, 1., 0.), std::invalid_argument);
  EXPECT_THROW(AqueductMPCAssemblyStorage(std::vector<double>{}), std::invalid_argument);
}
