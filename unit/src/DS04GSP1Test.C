#include "ADUtils.h"
#include "DS04GSP1.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <fstream>
#include <sstream>
namespace
{
using Nodes = DS04GSP1::Nodes;
using K8 = std::array<std::array<double, 8>, 8>;
constexpr double E = 2.97915e10, nu = .2, mu = E / (2 * (1 + nu)),
                 lam = E * nu / ((1 + nu) * (1 - 2 * nu));
Nodes nodes(bool distorted)
{
  if (!distorted)
    return {{{{0, 0}}, {{1, 0}}, {{1, 1}}, {{0, 1}}}};
  return {{{{41.0744375, 30.2808984}},
           {{40.5477539, 31.1235957}},
           {{39.777457, 30.5690957}},
           {{40.3079492, 29.7433555}}}};
}
K8 native(bool distorted)
{
  K8 K{};
  std::ifstream in(std::string("unit/data/ds04g-sp1/") + (distorted ? "e833" : "square") + ".mtx");
  if (!in)
    throw std::runtime_error("Missing frozen native matrix");
  std::string line;
  unsigned entries = 0;
  while (std::getline(in, line))
  {
    if (line.empty() || line[0] == '*')
      continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    unsigned n, d, m, e;
    double v;
    row >> n >> d >> m >> e >> v;
    if (!row || n < 1 || n > 4 || m < 1 || m > 4 || d < 1 || d > 2 || e < 1 || e > 2)
      throw std::runtime_error("Malformed native fixture");
    const auto i = 2 * (n - 1) + d - 1, j = 2 * (m - 1) + e - 1;
    K[i][j] = K[j][i] = 1e6 * v;
    ++entries;
  }
  if (entries != 36)
    throw std::runtime_error("Native coverage differs");
  return K;
}
K8 physical(const DS04GSP1::Geometry &g)
{
  K8 K{};
  for (unsigned j = 0; j < 8; ++j)
  {
    const unsigned n = j / 2, c = j % 2;
    const double ex = c == 0 ? g.gradient[n][0] : 0, ey = c == 1 ? g.gradient[n][1] : 0;
    const double xy = c == 0 ? g.gradient[n][1] : g.gradient[n][0];
    const double sx = lam * (ex + ey) + 2 * mu * ex, sy = lam * (ex + ey) + 2 * mu * ey,
                 tau = mu * xy;
    for (unsigned i = 0; i < 4; ++i)
    {
      K[2 * i][j] = g.area * (sx * g.gradient[i][0] + tau * g.gradient[i][1]);
      K[2 * i + 1][j] = g.area * (tau * g.gradient[i][0] + sy * g.gradient[i][1]);
    }
  }
  return K;
}
} // namespace
TEST(DS04GSP1, BothRealNativeMatricesAndLowStiffnessHGMatch)
{
  for (bool distorted : {false, true})
  {
    const auto x = nodes(distorted);
    const auto g = DS04GSP1::geometry(x);
    const auto N = native(distorted);
    const auto P = physical(g);
    double error = 0, norm = 0, herror = 0, hnorm = 0;
    for (unsigned i = 0; i < 8; ++i)
      for (unsigned j = 0; j < 8; ++j)
      {
        const double hg =
            (i % 2 == j % 2) ? mu * g.stiffness_per_shear * g.gamma[i / 2] * g.gamma[j / 2] : 0.;
        error += std::pow(P[i][j] + hg - N[i][j], 2);
        norm += N[i][j] * N[i][j];
        herror += std::pow(hg - (N[i][j] - P[i][j]), 2);
        hnorm += hg * hg;
      }
    EXPECT_LT(std::sqrt(error / norm), 1e-8);
    EXPECT_LT(std::sqrt(herror / hnorm), 1e-8);
  }
}
TEST(DS04GSP1, AffineAndRigidModesAreNullWithoutSquareNormalizationAssumption)
{
  for (bool distorted : {false, true})
  {
    const auto x = nodes(distorted);
    const auto g = DS04GSP1::geometry(x);
    double norm = 0;
    for (double z : g.gamma)
      norm += z * z;
    if (distorted)
      EXPECT_GT(std::abs(std::sqrt(norm) - 2), 1e-5);
    for (unsigned c = 0; c < 2; ++c)
    {
      std::array<double, 4> u;
      for (unsigned j = 0; j < 4; ++j)
        u[j] = .3 + .002 * x[j][0] - .004 * x[j][1];
      for (unsigned i = 0; i < 4; ++i)
        EXPECT_NEAR(DS04GSP1::hourglass(g, u, i, mu), 0, 1e-7);
    }
  }
}
TEST(DS04GSP1, HGFiniteRigidSuperpositionPreservesEnergyAndRotatesForce)
{
  const auto x = nodes(true);
  const auto g = DS04GSP1::geometry(x);
  std::array<double, 4> u{{.001, -.003, .002, .005}}, v{{-.002, .003, .001, -.001}};
  const double a = .7, c = std::cos(a), s = std::sin(a);
  std::array<double, 4> ru, rv;
  for (unsigned i = 0; i < 4; ++i)
  {
    ru[i] = c * (x[i][0] + u[i]) - s * (x[i][1] + v[i]) + 4 - x[i][0];
    rv[i] = s * (x[i][0] + u[i]) + c * (x[i][1] + v[i]) - 3 - x[i][1];
  }
  for (unsigned i = 0; i < 4; ++i)
  {
    const auto fx = DS04GSP1::hourglass(g, u, i, mu), fy = DS04GSP1::hourglass(g, v, i, mu);
    EXPECT_NEAR(DS04GSP1::hourglass(g, ru, i, mu), c * fx - s * fy, 1e-5);
    EXPECT_NEAR(DS04GSP1::hourglass(g, rv, i, mu), s * fx + c * fy, 1e-5);
  }
}
TEST(DS04GSP1, MassCompensationRestoresSquareMatrixAndDistortedRowsAreNotZero)
{
  const auto g = DS04GSP1::geometry(nodes(false));
  constexpr double expected[4][4] = {{4, 2, 1, 2}, {2, 4, 2, 1}, {1, 2, 4, 2}, {2, 1, 2, 4}};
  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 4; ++j)
      EXPECT_NEAR(g.mass_correction[i][j] + g.area / 16, expected[i][j] / 36., 1e-14);
  const auto d = DS04GSP1::geometry(nodes(true));
  double sum = 0, maxrow = 0;
  for (const auto &row : d.mass_correction)
  {
    double s = 0;
    for (double z : row)
      s += z;
    sum += s;
    maxrow = std::max(maxrow, std::abs(s));
  }
  EXPECT_NEAR(sum, 0, 1e-14);
  EXPECT_GT(maxrow, 1e-5);
}
TEST(DS04GSP1, CurrentGeometryGravityADAndFiniteDifferenceAgree)
{
  auto x = nodes(true);
  x[0][0] += .03;
  x[2][1] -= .02;
  std::array<std::array<ADReal, 2>, 4> ad;
  for (unsigned i = 0; i < 4; ++i)
    for (unsigned c = 0; c < 2; ++c)
    {
      ad[i][c] = x[i][c];
      Moose::derivInsert(ad[i][c].derivatives(), 2 * i + c, 1.);
    }
  const auto r = DS04GSP1::gravityCorrection(ad);
  for (unsigned k = 0; k < 8; ++k)
  {
    auto plus = x, minus = x;
    plus[k / 2][k % 2] += 1e-6;
    minus[k / 2][k % 2] -= 1e-6;
    const auto p = DS04GSP1::gravityCorrection(plus), m = DS04GSP1::gravityCorrection(minus);
    for (unsigned i = 0; i < 4; ++i)
      EXPECT_NEAR(r[i].derivatives()[k], (p[i] - m[i]) / 2e-6, 2e-8);
  }
  const auto delta = DS04GSP1::gravityCorrection(x);
  double s = 0;
  for (auto z : delta)
    s += z;
  EXPECT_NEAR(s, 0, 1e-13);
}
TEST(DS04GSP1, HGADCurrentTangentAndRetryHistoryDoNotMutateOld)
{
  const auto g = DS04GSP1::geometry(nodes(true));
  std::array<ADReal, 4> u;
  for (unsigned i = 0; i < 4; ++i)
  {
    u[i] = .003 * i;
    Moose::derivInsert(u[i].derivatives(), i, 1.);
  }
  for (double h : {.002, .001, .00025})
  {
    const auto w = DS04GSP1::history(-.05, .003604, h);
    for (unsigned i = 0; i < 4; ++i)
    {
      const auto f = w.current * DS04GSP1::hourglass(g, u, i, mu);
      for (unsigned j = 0; j < 4; ++j)
        EXPECT_NEAR(f.derivatives()[j],
                    w.current * mu * g.stiffness_per_shear * g.gamma[i] * g.gamma[j], 1e-4);
    }
    EXPECT_NEAR(w.current + w.old + w.older, 1., 1e-14);
  }
  EXPECT_EQ(DS04GSP1::history(-.05, .003604, .002).current,
            DS04GSP1::history(-.05, .003604, .002).current);
}
TEST(DS04GSP1, NewmarkAccelerationUsesTrialIntervalAndAcceptedHistory)
{
  const double u = .1, v = .03, a = -.7;
  for (double h : {.002, .001, .00025})
    EXPECT_NEAR(DS04GSP1::trialAcceleration(u + h * v + .5 * h * h * a, u, v, a, h, .275625), a,
                1e-8);
  EXPECT_THROW(DS04GSP1::trialAcceleration(u, u, v, a, 0., .275625), std::invalid_argument);
}
TEST(DS04GSP1, InvalidAndTranslatedGeometry)
{
  auto x = nodes(true);
  const auto a = DS04GSP1::geometry(x);
  for (auto &p : x)
  {
    p[0] += 10000;
    p[1] -= 10000;
  }
  const auto b = DS04GSP1::geometry(x);
  EXPECT_NEAR(a.stiffness_per_shear, b.stiffness_per_shear, 1e-12);
  std::swap(x[0], x[1]);
  EXPECT_THROW(DS04GSP1::geometry(x), std::invalid_argument);
}
