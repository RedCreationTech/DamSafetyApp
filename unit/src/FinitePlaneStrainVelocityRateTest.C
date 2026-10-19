#include "FinitePlaneStrainVelocityRate.h"
#include "gtest/gtest.h"
#include "ADUtils.h"
#include <algorithm>
#include <iostream>

namespace
{
constexpr Real kYoungs = 2.97915e10, kNu = .2, kZeta = .003604;
constexpr Real kBeta = .275625, kGamma = .55, kAlpha = -.05, kDt = .002;
RankTwoTensor
tensor(Real xx, Real xy, Real yx, Real yy)
{
  RankTwoTensor a;
  a(0, 0) = xx;
  a(0, 1) = xy;
  a(1, 0) = yx;
  a(1, 1) = yy;
  return a;
}
}

TEST(FinitePlaneStrainVelocityRate, ReusesNewmarkAcceptedHistoryAcrossRetries)
{
  const Real u0 = .02, v0 = .07, a0 = -.3;
  for (Real h : {.01, .005, .0037})
  {
    const Real u = u0 + h * v0 + .5 * h * h * a0;
    EXPECT_NEAR(
        ElasticVelocityRate::trialVelocity(u, u0, v0, a0, h, kBeta, kGamma), v0 + h * a0, 1e-13);
  }
  const auto before = ElasticVelocityRate::trialVelocity(.031, u0, v0, a0, .002, kBeta, kGamma);
  (void)ElasticVelocityRate::trialVelocity(.035, u0, v0, a0, .004, kBeta, kGamma);
  EXPECT_DOUBLE_EQ(before,
                   ElasticVelocityRate::trialVelocity(.031, u0, v0, a0, .002, kBeta, kGamma));
}

TEST(FinitePlaneStrainVelocityRate, PlaneStrainRetainsNonzeroOutOfPlaneStress)
{
  const auto f = FinitePlaneStrainVelocityRate::deformation(RankTwoTensor());
  const auto rate = tensor(.001, .0002, .0002, -.0003);
  const auto stress = FinitePlaneStrainVelocityRate::stress(rate, f, kYoungs, kNu);
  const Real lambda = kYoungs * kNu / ((1 + kNu) * (1 - 2 * kNu)), mu = kYoungs / (2 * (1 + kNu));
  EXPECT_NEAR(stress(2, 2), lambda * .0007, 1e-8);
  EXPECT_NEAR(stress(0, 0), 2 * mu * .001 + lambda * .0007, 1e-8);
  EXPECT_NEAR(stress(0, 1), 2 * mu * .0002, 1e-8);
}

TEST(FinitePlaneStrainVelocityRate, InstantRigidRotationHasZeroDampingStress)
{
  for (Real angle : {0., .3, 1.2})
  {
    auto f = tensor(std::cos(angle), -std::sin(angle), std::sin(angle), std::cos(angle));
    f(2, 2) = 1.;
    const auto omega = tensor(0., -.7, .7, 0.);
    const auto stress = FinitePlaneStrainVelocityRate::stress(omega * f, f, kYoungs, kNu);
    EXPECT_NEAR(stress.L2norm(), 0., 2e-5);
  }
}

TEST(FinitePlaneStrainVelocityRate, StaticPredeformationWithZeroVelocityIsAtRest)
{
  const auto g = tensor(.12, .07, -.02, -.03);
  const auto p = FinitePlaneStrainVelocityRate::nominalStress(
      g, g, RankTwoTensor(), RankTwoTensor(), kDt, kBeta, kGamma, kAlpha, kZeta, kYoungs, kNu);
  EXPECT_DOUBLE_EQ(p.L2norm(), 0.);
}

TEST(FinitePlaneStrainVelocityRate, HHTUsesAcceptedOldSpatialRateAndCurrentPullback)
{
  const auto g = tensor(.12, .03, -.02, -.04), old = tensor(.1, .02, -.01, -.03);
  const auto v = tensor(.004, .002, -.003, .001), a = tensor(.01, .02, .03, -.02);
  auto f = g;
  f.addIa(1.);
  auto old_f = old;
  old_f.addIa(1.);
  const auto acceleration = ((g - old) / (kDt * kDt) - v / kDt - a * (.5 - kBeta)) / kBeta;
  const auto velocity = v + a * (kDt * (1 - kGamma)) + acceleration * (kDt * kGamma);
  const auto l = velocity * f.inverse(), lo = v * old_f.inverse();
  const auto d = (l + l.transpose()) * .5, d_old = (lo + lo.transpose()) * .5;
  const Real mu = kYoungs / (2 * (1 + kNu)), lambda = kYoungs * kNu / ((1 + kNu) * (1 - 2 * kNu));
  auto s = d * (2 * mu);
  s.addIa(lambda * d.trace());
  auto so = d_old * (2 * mu);
  so.addIa(lambda * d_old.trace());
  const auto expected =
      (((s * (1 + kAlpha) - so * kAlpha) * kZeta) * f.inverse().transpose()) * f.det();
  const auto actual = FinitePlaneStrainVelocityRate::nominalStress(
      g, old, v, a, kDt, kBeta, kGamma, kAlpha, kZeta, kYoungs, kNu);
  EXPECT_LE((actual - expected).L2norm(), 1e-12 * expected.L2norm());
}

TEST(FinitePlaneStrainVelocityRate, FullFinitePullbackADMatchesFourGradientPerturbations)
{
  const auto g = tensor(.12, .03, -.02, -.04), old = tensor(.119, .029, -.019, -.039);
  const auto v = tensor(.004, .002, -.003, .001), a = tensor(.01, .02, .03, -.02);
  ADRankTwoTensor ad_g = g;
  for (unsigned int k = 0; k < 4; ++k)
    Moose::derivInsert(ad_g(k / 2, k % 2).derivatives(), k, 1.);
  const auto ad_p = FinitePlaneStrainVelocityRate::nominalStress(
      ad_g, old, v, a, kDt, kBeta, kGamma, kAlpha, kZeta, kYoungs, kNu);
  for (unsigned int k = 0; k < 4; ++k)
  {
    auto plus = g, minus = g;
    constexpr Real h = 1e-7;
    plus(k / 2, k % 2) += h;
    minus(k / 2, k % 2) -= h;
    const auto p = FinitePlaneStrainVelocityRate::nominalStress(
        plus, old, v, a, kDt, kBeta, kGamma, kAlpha, kZeta, kYoungs, kNu);
    const auto m = FinitePlaneStrainVelocityRate::nominalStress(
        minus, old, v, a, kDt, kBeta, kGamma, kAlpha, kZeta, kYoungs, kNu);
    for (unsigned int i = 0; i < 2; ++i)
      for (unsigned int j = 0; j < 2; ++j)
      {
        const Real fd = (p(i, j) - m(i, j)) / (2 * h);
        EXPECT_NEAR(ad_p(i, j).derivatives()[k], fd, 2e-7 * std::max(1., std::abs(fd)));
      }
  }
}

TEST(FinitePlaneStrainVelocityRate, ZeroDampingIsExactlyZeroIncludingDerivative)
{
  ADRankTwoTensor g = tensor(.12, .03, -.02, -.04);
  Moose::derivInsert(g(0, 0).derivatives(), 0, 1.);
  const auto p = FinitePlaneStrainVelocityRate::nominalStress(g,
                                                              RankTwoTensor(),
                                                              RankTwoTensor(),
                                                              RankTwoTensor(),
                                                              kDt,
                                                              kBeta,
                                                              kGamma,
                                                              kAlpha,
                                                              0.,
                                                              kYoungs,
                                                              kNu);
  EXPECT_DOUBLE_EQ(MetaPhysicL::raw_value(p.L2norm()), 0.);
  EXPECT_DOUBLE_EQ(p(0, 0).derivatives()[0], 0.);
}

TEST(FinitePlaneStrainVelocityRate, RejectsInvertedCurrentAndOldConfigurations)
{
  const auto bad = tensor(-2., 0., 0., 0.), good = RankTwoTensor();
  EXPECT_THROW(FinitePlaneStrainVelocityRate::nominalStress(
                   bad, good, good, good, kDt, kBeta, kGamma, kAlpha, kZeta, kYoungs, kNu),
               std::invalid_argument);
  EXPECT_THROW(FinitePlaneStrainVelocityRate::nominalStress(
                   good, bad, good, good, kDt, kBeta, kGamma, kAlpha, kZeta, kYoungs, kNu),
               std::invalid_argument);
}

namespace
{
using Nodes = std::array<std::array<Real, 2>, 4>;
template <typename T>
std::array<T, 8>
assembledForce(const Nodes & nodes,
               const std::array<T, 8> & u,
               const std::array<Real, 8> & old,
               const std::array<Real, 8> & v,
               const std::array<Real, 8> & a,
               Real h)
{
  std::array<T, 8> force{};
  const Real g = 1 / std::sqrt(3.0);
  const std::array<Real, 4> sx{{-1, 1, 1, -1}}, sy{{-1, -1, 1, 1}};
  for (Real xi : {-g, g})
    for (Real eta : {-g, g})
    {
      RankTwoTensor map;
      map(2, 2) = 1;
      std::array<std::array<Real, 2>, 4> b;
      for (unsigned int n = 0; n < 4; ++n)
      {
        b[n] = {{sx[n] * (1 + sy[n] * eta) / 4, sy[n] * (1 + sx[n] * xi) / 4}};
        for (unsigned int i = 0; i < 2; ++i)
          for (unsigned int j = 0; j < 2; ++j)
            map(i, j) += nodes[n][i] * b[n][j];
      }
      const auto inverse = map.inverse();
      for (auto & row : b)
      {
        const auto original = row;
        for (unsigned int j = 0; j < 2; ++j)
          row[j] = original[0] * inverse(0, j) + original[1] * inverse(1, j);
      }
      RankTwoTensorTempl<T> gradient;
      RankTwoTensor old_gradient, velocity, acceleration;
      for (unsigned int n = 0; n < 4; ++n)
        for (unsigned int i = 0; i < 2; ++i)
          for (unsigned int j = 0; j < 2; ++j)
          {
            gradient(i, j) += u[2 * n + i] * b[n][j];
            old_gradient(i, j) += old[2 * n + i] * b[n][j];
            velocity(i, j) += v[2 * n + i] * b[n][j];
            acceleration(i, j) += a[2 * n + i] * b[n][j];
          }
      const auto p = FinitePlaneStrainVelocityRate::nominalStress(gradient,
                                                                  old_gradient,
                                                                  velocity,
                                                                  acceleration,
                                                                  h,
                                                                  kBeta,
                                                                  kGamma,
                                                                  kAlpha,
                                                                  kZeta,
                                                                  kYoungs,
                                                                  kNu);
      for (unsigned int n = 0; n < 4; ++n)
        for (unsigned int i = 0; i < 2; ++i)
          for (unsigned int j = 0; j < 2; ++j)
            force[2 * n + i] += p(i, j) * b[n][j] * map.det();
    }
  return force;
}
}

TEST(FinitePlaneStrainVelocityRate, FourGaussEightDofAssemblyMatchesFiniteDifference)
{
  const std::array<Nodes, 3> geometries{
      {{{{{0, 0}}, {{2, 0}}, {{2, 1}}, {{0, 1}}}},
       {{{{.1, -.2}}, {{2.2, .1}}, {{1.8, 1.7}}, {{-.3, 1.2}}}},
       // Original DS04G CPE4R 1777 corner coordinates, mm converted to m.
       {{{{10, 80}}, {{9, 80}}, {{9.29290527, 78.6176484}}, {{10.5266855, 79.1573047}}}}}};
  Real maximum = 0;
  for (const auto & nodes : geometries)
    for (Real step : {.002, .0013})
      for (Real stretch : {0., .01, .08, .18})
      {
        std::array<Real, 8> u, old, v, a;
        std::array<ADReal, 8> ad_u;
        for (unsigned int k = 0; k < 8; ++k)
        {
          old[k] = stretch * nodes[k / 2][k % 2] + .0003 * (k % 3);
          v[k] = .004 * (static_cast<int>(k) % 3 - 1);
          a[k] = .02 * (static_cast<int>(k) % 4 - 2);
          u[k] = old[k] + step * v[k] + .5 * step * step * a[k] + .00001 * (k % 2);
          ad_u[k] = u[k];
          Moose::derivInsert(ad_u[k].derivatives(), k, 1.);
        }
        const auto ad = assembledForce(nodes, ad_u, old, v, a, step);
        Real scale = 1, error = 0;
        for (unsigned int k = 0; k < 8; ++k)
        {
          auto plus = u, minus = u;
          const Real perturbation = 2e-6 * std::max(1., std::abs(u[k]));
          plus[k] += perturbation;
          minus[k] -= perturbation;
          const auto fp = assembledForce(nodes, plus, old, v, a, step);
          const auto fm = assembledForce(nodes, minus, old, v, a, step);
          for (unsigned int i = 0; i < 8; ++i)
          {
            const Real tangent = ad[i].derivatives()[k];
            const Real fd = (fp[i] - fm[i]) / (2 * perturbation);
            scale = std::max(scale, std::abs(tangent));
            error = std::max(error, std::abs(tangent - fd));
          }
        }
        maximum = std::max(maximum, error / scale);
        EXPECT_LT(error / scale, 1e-7);
      }
  std::cout << "24 four-Gauss 8DOF formula assemblies max relative AD/FD error=" << maximum << '\n';
}
