#include "AbaqusCDPSubstepIntegrator.h"
#include "CombinedDamageViscosity.h"
#include "MooseTypes.h"
#include "ADUtils.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
using Core = AbaqusCDPStateIntegrator;
using Tensor = Core::SymmetricTensor;
constexpr double E = 3.04e10, nu = .2, mu = .0005;
CDPMaterialTable table()
{
  const std::string p = "test/tests/cdp_material_table/data/";
  return {p + "compression_hardening.csv", p + "compression_damage.csv",
          p + "tension_stiffening.csv", p + "tension_damage.csv", E};
}
AbaqusCDPLocalIntegrator::Parameters localParameters(unsigned int iterations = 40)
{
  return {E, nu, 36.31, .1, 1.16, .667, iterations, 1e-9, 1e-7, 1e-6};
}
Core::Parameters parameters(bool combined, double viscosity = mu, double wt = 1., double wc = 1.)
{
  return {wt, wc, viscosity, 1e-12, combined};
}
Tensor strain(double stress)
{
  return {stress / E, -nu * stress / E, -nu * stress / E, 0., 0., 0.};
}
void sameBranches(const Core::Result & a, const Core::Result & b)
{
  EXPECT_EQ(a.state.backbone.plastic_strain, b.state.backbone.plastic_strain);
  EXPECT_DOUBLE_EQ(a.state.backbone.tensile_equivalent_plastic_strain,
                   b.state.backbone.tensile_equivalent_plastic_strain);
  EXPECT_DOUBLE_EQ(a.state.backbone.compressive_equivalent_plastic_strain,
                   b.state.backbone.compressive_equivalent_plastic_strain);
  EXPECT_EQ(a.state.viscous_plastic_strain, b.state.viscous_plastic_strain);
  EXPECT_DOUBLE_EQ(a.state.viscous_tension_damage, b.state.viscous_tension_damage);
  EXPECT_DOUBLE_EQ(a.state.viscous_compression_damage, b.state.viscous_compression_damage);
  EXPECT_EQ(a.viscous_effective_stress, b.viscous_effective_stress);
  EXPECT_DOUBLE_EQ(a.damage.tension_weight, b.damage.tension_weight);
}
Core::State damagedState(const CDPMaterialTable & t)
{
  Core::State s;
  s.backbone.tensile_equivalent_plastic_strain = 4.1e-5;
  s.backbone.compressive_equivalent_plastic_strain = 7.1e-5;
  s.viscous_tension_damage = .5 * t.responseByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::TENSION, s.backbone.tensile_equivalent_plastic_strain).damage.value;
  s.viscous_compression_damage = .5 * t.responseByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::COMPRESSION, s.backbone.compressive_equivalent_plastic_strain).damage.value;
  s.viscous_combined_damage = .12;
  return s;
}
ADReal seeded(double value, double derivative, unsigned int index = 0)
{
  ADReal x = value;
  Moose::derivInsert(x.derivatives(), index, derivative);
  return x;
}
}

TEST(CombinedDamageViscosity, ActualADIncludesTargetAndOldHistoryAndRecovers)
{
  for (const double h : {0., 2.5e-4, .01})
    for (const double viscosity : {0., mu})
    {
      ADReal old = seeded(.4, 1., 0), target = seeded(.1, 1., 1);
      const ADReal d = CombinedDamageViscosity::update(old, target, h, viscosity);
      const double a = viscosity == 0. ? 1. : h / (viscosity + h);
      EXPECT_NEAR(d.value(), (1. - a) * .4 + a * .1, 1e-16);
      EXPECT_NEAR(d.derivatives()[0], 1. - a, 1e-15);
      EXPECT_NEAR(d.derivatives()[1], a, 1e-15);
    }
}

TEST(CombinedDamageViscosity, DefaultLegacyIsIndependentOfAppendedHistory)
{
  const auto t = table();
  const AbaqusCDPLocalIntegrator local(t, localParameters());
  const Core implicit(local, {1., 1., mu, 1e-12});
  const Core explicit_legacy(local, parameters(false));
  auto old = damagedState(t);
  const auto target = strain(1e6);
  const auto a = implicit.integrateLinearized(target, .01, old);
  old.viscous_combined_damage = .77;
  const auto b = explicit_legacy.integrateLinearized(target, .01, old);
  EXPECT_EQ(a.result.cauchy_stress, b.result.cauchy_stress);
  EXPECT_EQ(a.derivative, b.derivative);
  EXPECT_DOUBLE_EQ(a.result.damage.damage, a.result.legacy_damage_reconstruction);
  EXPECT_DOUBLE_EQ(a.result.state.viscous_combined_damage, a.result.damage.damage);
  for (unsigned int i = 0; i < 6; ++i) EXPECT_EQ(a.derivative[22][i], 0.);
}

TEST(CombinedDamageViscosity, FixedPathLeavesBackbonePlasticAndBranchesUnchanged)
{
  const auto t = table();
  const AbaqusCDPLocalIntegrator local(t, localParameters());
  const Core legacy(local, parameters(false)), candidate(local, parameters(true));
  Core::State a, b;
  const double ft = t.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.).stress.value;
  const double fc = t.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::COMPRESSION, 0.).stress.value;
  for (const double s : {.5 * ft, 1.05 * ft, .4 * ft, -.5 * fc, -1.05 * fc, -.2 * fc, .3 * ft})
  {
    const auto ra = legacy.integrate(strain(s), .00025, a);
    const auto rb = candidate.integrate(strain(s), .00025, b);
    sameBranches(ra, rb);
    a = ra.state;
    b = rb.state;
  }
}

TEST(CombinedDamageViscosity, ZeroViscosityAndZeroTimeBoundaries)
{
  const auto t = table();
  const AbaqusCDPLocalIntegrator local(t, localParameters());
  const Core legacy(local, parameters(false, 0.)), candidate(local, parameters(true, 0.));
  const auto old = damagedState(t);
  for (const double h : {0., .01})
  {
    const auto a = legacy.integrateLinearized(strain(1e6), h, old);
    const auto b = candidate.integrateLinearized(strain(1e6), h, old);
    sameBranches(a.result, b.result);
    for (unsigned int i = 0; i < 6; ++i)
    {
      EXPECT_NEAR(a.result.cauchy_stress[i], b.result.cauchy_stress[i], 1e-8);
      for (unsigned int j = 0; j < 23; ++j)
        EXPECT_NEAR(a.derivative[j][i], b.derivative[j][i], 1e-5);
    }
  }
  const Core viscous(local, parameters(true));
  const auto hold = viscous.integrateLinearized(strain(-1e6), 0., old);
  EXPECT_DOUBLE_EQ(hold.result.state.viscous_combined_damage, old.viscous_combined_damage);
  EXPECT_EQ(hold.derivative[22][22], 1.);
}

TEST(CombinedDamageViscosity, RecoveryDecreasesCombinedScalarWithoutErasingBranchDamage)
{
  const auto t = table();
  const AbaqusCDPLocalIntegrator local(t, localParameters());
  const Core candidate(local, parameters(true));
  auto old = damagedState(t);
  old.viscous_combined_damage = .8;
  const auto result = candidate.integrate(strain(-1e6), .001, old);
  EXPECT_LT(result.state.viscous_combined_damage, old.viscous_combined_damage);
  EXPECT_GE(result.state.viscous_tension_damage, old.viscous_tension_damage);
  EXPECT_GE(result.state.viscous_compression_damage, old.viscous_compression_damage);
  EXPECT_NEAR(result.damage.stiffness_factor, 1. - result.state.viscous_combined_damage, 1e-16);
  old.viscous_combined_damage = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(candidate.integrate(strain(-1e6), .001, old), std::runtime_error);
}

TEST(CombinedDamageViscosity, ActualADChecksFullStateColumnsIncludingWeightAndHistory)
{
  const auto t = table();
  const AbaqusCDPLocalIntegrator local(t, localParameters());
  // Mixed diagonal stress: distinct principal values, no eigenvalue crossing.
  // The unchanged backbone Jacobian supplies seeds; the new combination,
  // scalar history, effective-stress weight and stress product use real AD.
  for (const auto recovery : {std::array<double, 2>{1., 1.}, {.63, .27}})
  {
    const Core candidate(local, parameters(true, mu, recovery[0], recovery[1]));
    const auto old = damagedState(t);
    const Tensor eps = {1e-5, -1.5e-5, 3e-6, 0., 0., 0.};
    const double h = .00025;
    const auto bone = local.integrateLinearized(eps, old.backbone);
    const auto result = candidate.integrateLinearized(eps, h, old);
    ASSERT_FALSE(bone.result.plastic);
    for (unsigned int col = 0; col < Core::transition_size; ++col)
    {
      std::array<ADReal, 6> epv, sig;
      std::array<ADReal, 8> bone_state;
      for (unsigned int k = 0; k < 8; ++k)
      {
        const double v = k < 6 ? bone.result.state.plastic_strain[k] :
            (k == 6 ? bone.result.state.tensile_equivalent_plastic_strain :
                      bone.result.state.compressive_equivalent_plastic_strain);
        bone_state[k] = seeded(v, col < 14 ? bone.derivative[col][6 + k] : 0.);
      }
      for (unsigned int k = 0; k < 6; ++k)
        epv[k] = CombinedDamageViscosity::update(
            seeded(old.viscous_plastic_strain[k], col == 14 + k ? 1. : 0.),
            bone_state[k], h, mu);
      const double G = E / (2. * (1. + nu));
      const double lambda = E * nu / ((1. + nu) * (1. - 2. * nu));
      std::array<ADReal, 6> elastic;
      for (unsigned int k = 0; k < 6; ++k)
        elastic[k] = seeded(eps[k], col == k ? 1. : 0.) - epv[k];
      for (unsigned int k = 0; k < 6; ++k)
        sig[k] = 2. * G * elastic[k] + (k < 3 ? lambda * (elastic[0] + elastic[1] + elastic[2]) : ADReal(0.));
      ASSERT_GT(sig[0].value(), 0.);
      ASSERT_LT(sig[1].value(), 0.);
      ASSERT_GT(sig[2].value(), 0.);
      const ADReal rv = (sig[0] + sig[2]) / (sig[0] - sig[1] + sig[2]);
      const auto bt = local.materialResponse(CDPMaterialTable::Branch::TENSION,
          bone.result.state.tensile_equivalent_plastic_strain);
      const auto bc = local.materialResponse(CDPMaterialTable::Branch::COMPRESSION,
          bone.result.state.compressive_equivalent_plastic_strain);
      ADReal dt = seeded(bt.damage.value, bt.damage.right_derivative * bone_state[6].derivatives()[0]);
      ADReal dc = seeded(bc.damage.value, bc.damage.right_derivative * bone_state[7].derivatives()[0]);
      const ADReal target = 1. - (1. - (1. - recovery[0] * rv) * dc) *
                                   (1. - (1. - recovery[1] * (1. - rv)) * dt);
      const ADReal dv = CombinedDamageViscosity::update(
          seeded(old.viscous_combined_damage, col == 22 ? 1. : 0.), target, h, mu);
      EXPECT_NEAR(dv.value(), result.result.state.viscous_combined_damage, 1e-14);
      EXPECT_NEAR(dv.derivatives()[0], result.derivative[col][22],
                  2e-6 * std::max(1., std::abs(dv.derivatives()[0])));
      for (unsigned int k = 0; k < 6; ++k)
      {
        const ADReal stress = (1. - dv) * sig[k];
        EXPECT_NEAR(stress.value(), result.result.cauchy_stress[k], 1e-7);
        EXPECT_NEAR(stress.derivatives()[0], result.derivative[col][k],
                    2e-6 * std::max(1., std::abs(stress.derivatives()[0])));
      }
    }
  }
}

TEST(CombinedDamageViscosity, PlasticSubstepTangentMatchesIndependentPerturbations)
{
  const auto t = table();
  const AbaqusCDPLocalIntegrator local(t, localParameters());
  const Core candidate(local, parameters(true));
  const double ft = t.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.).stress.value;
  const double fc = t.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::COMPRESSION, 0.).stress.value;
  for (const double stress : {1.05 * ft, -1.05 * fc})
  {
    auto target = strain(stress);
    // Distinct transverse strains avoid the repeated plastic eigenvalues of
    // exact uniaxial loading, where a centered derivative is not unique.
    target[1] *= .7;
    target[2] *= 1.3;
    target[3] = .02 * target[0];
    const AbaqusCDPSubstepIntegrator steps(candidate, {16, std::abs(target[0]) / 3., 1e-8});
    const auto linear = steps.integrateLinearized({}, target, .001, {});
    ASSERT_EQ(linear.result.accepted_substeps, 4u);
    for (const Tensor direction : {Tensor{1., -.2, -.2, .1, 0., 0.}, Tensor{.2, .3, -.4, 0., .05, .03}})
    {
      SCOPED_TRACE(std::to_string(stress) + ":" + std::to_string(direction[0]));
      const auto predicted = steps.applyTangent(linear.algorithmic_tangent, direction);
      const auto measured = steps.directionalDerivative({}, target, .001, {}, direction, 1e-9);
      for (unsigned int k = 0; k < 6; ++k)
        EXPECT_NEAR(predicted[k], measured[k], .01 * std::max(1., std::abs(measured[k])));
    }
  }
}

TEST(CombinedDamageViscosity, FailedPartialPartitionsResetCombinedHistoryAndRetryIsDeterministic)
{
  const auto t = table();
  const AbaqusCDPLocalIntegrator failing(t, localParameters(0)), healthy(t, localParameters());
  const Core bad(failing, parameters(true)), good(healthy, parameters(true));
  const AbaqusCDPSubstepIntegrator steps(bad, {8, 0., 1e-8});
  Core::State old;
  old.viscous_combined_damage = .2;
  const auto snapshot = old;
  const double ft = t.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.).stress.value;
  AbaqusCDPSubstepIntegrator::Trace trace;
  EXPECT_THROW(steps.integrateLinearized({}, strain(2.1 * ft), .001, old, &trace), std::runtime_error);
  bool partial_success = false;
  for (const auto & e : trace)
  {
    if (e.substep == 1) EXPECT_DOUBLE_EQ(e.old_state.viscous_combined_damage, .2);
    if (e.succeeded)
    {
      partial_success = true;
      EXPECT_LT(e.result.state.viscous_combined_damage, .2);
    }
  }
  EXPECT_TRUE(partial_success);
  EXPECT_DOUBLE_EQ(old.viscous_combined_damage, snapshot.viscous_combined_damage);
  const auto clean = good.integrateLinearized(strain(.5 * ft), .0005, snapshot);
  const auto retry = good.integrateLinearized(strain(.5 * ft), .0005, old);
  EXPECT_EQ(clean.derivative, retry.derivative);
  EXPECT_EQ(clean.result.cauchy_stress, retry.result.cauchy_stress);
  EXPECT_DOUBLE_EQ(clean.result.state.viscous_combined_damage, retry.result.state.viscous_combined_damage);
}
