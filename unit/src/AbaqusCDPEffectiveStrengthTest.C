#include "AbaqusCDPStateIntegrator.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace
{
constexpr double modulus = 2.97915e10;
constexpr double poisson = 0.2;
const std::string fixture = "test/tests/cdp_material_table/uniaxial_tension_20260830/";
using Branch = CDPMaterialTable::Branch;
using Tensor = AbaqusCDPFormula::SymmetricTensor;

CDPMaterialTable
expertTable()
{
  return CDPMaterialTable(fixture + "compression_hardening.csv",
                         fixture + "compression_damage.csv",
                         fixture + "tension_stiffening.csv",
                         fixture + "tension_damage.csv",
                         modulus);
}

AbaqusCDPLocalIntegrator::Parameters
localParameters(const bool automatic = true)
{
  auto result = AbaqusCDPLocalIntegrator::Parameters{modulus, poisson, 36.0, 0.1, 1.16, 0.667};
  result.use_automatic_differentiation_jacobian = automatic;
  return result;
}

AbaqusCDPStateIntegrator::State
damagedState(const CDPMaterialTable & table, const Branch branch, const double kappa)
{
  AbaqusCDPStateIntegrator::State state;
  const auto response = table.responseByEquivalentPlasticStrain(branch, kappa);
  if (branch == Branch::TENSION)
  {
    state.backbone.tensile_equivalent_plastic_strain = kappa;
    state.backbone.plastic_strain[0] = kappa;
    state.viscous_tension_damage = response.damage.value;
  }
  else
  {
    state.backbone.compressive_equivalent_plastic_strain = kappa;
    state.backbone.plastic_strain[0] = -kappa;
    state.viscous_compression_damage = response.damage.value;
  }
  state.viscous_plastic_strain = state.backbone.plastic_strain;
  return state;
}

Tensor
strainAtEffectiveStress(const AbaqusCDPStateIntegrator::State & state, const double stress)
{
  auto strain = state.backbone.plastic_strain;
  strain[0] += stress / modulus;
  strain[1] -= poisson * stress / modulus;
  strain[2] -= poisson * stress / modulus;
  return strain;
}
}

// These are admissible frozen-history states, not a substitute for integrating
// an entire tensile/compressive loading path or for an Abaqus comparison.
TEST(AbaqusCDPEffectiveStrength, DamagedTableStatesStayOnEffectiveYieldSurface)
{
  const auto table = expertTable();
  const AbaqusCDPLocalIntegrator local(table, localParameters());
  for (const auto branch : {Branch::TENSION, Branch::COMPRESSION})
  {
    const auto & kappas = table.equivalentPlasticStrain(branch);
    for (std::size_t row = 1; row < kappas.size(); ++row)
    {
      SCOPED_TRACE(std::to_string(static_cast<int>(branch)) + ":" + std::to_string(row));
      const double nominal = table.stressValues(branch)[row];
      const double damage = table.response(branch, table.stressAbscissa(branch)[row]).damage.value;
      const double sign = branch == Branch::TENSION ? 1.0 : -1.0;
      const double effective = sign * nominal / (1.0 - damage);
      const auto old = damagedState(table, branch, kappas[row]);
      const auto result = local.integrate(strainAtEffectiveStress(old, effective), old.backbone);
      EXPECT_FALSE(result.plastic);
      EXPECT_EQ(result.state.plastic_strain, old.backbone.plastic_strain);
      EXPECT_NEAR(result.trial_yield / std::abs(effective), 0.0, 1.0e-10);
      EXPECT_NEAR(result.effective_stress[0], effective, std::abs(effective) * 1.0e-12);
    }
  }
}

TEST(AbaqusCDPEffectiveStrength, DamagedElasticUnloadingDoesNotTriggerPlasticCorrection)
{
  const auto table = expertTable();
  const AbaqusCDPLocalIntegrator local(table, localParameters());
  for (const auto branch : {Branch::TENSION, Branch::COMPRESSION})
  {
    const double kappa = table.equivalentPlasticStrain(branch).back();
    const auto sample = table.responseByEquivalentPlasticStrain(branch, kappa);
    const auto old = damagedState(table, branch, kappa);
    const double sign = branch == Branch::TENSION ? 1.0 : -1.0;
    const double stress = sign * 0.9 * sample.stress.value / (1.0 - sample.damage.value);
    const auto result = local.integrate(strainAtEffectiveStress(old, stress), old.backbone);
    EXPECT_FALSE(result.plastic);
    EXPECT_LT(result.trial_yield, 0.0);
    EXPECT_EQ(result.state.plastic_strain, old.backbone.plastic_strain);
  }
}

TEST(AbaqusCDPEffectiveStrength, ZeroViscosityRecoversNominalTableStressExactlyOnce)
{
  const auto table = expertTable();
  const AbaqusCDPLocalIntegrator local(table, localParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, {0.0, 1.0, 0.0, 1.0e-12});
  for (const auto branch : {Branch::TENSION, Branch::COMPRESSION})
  {
    const auto & kappas = table.equivalentPlasticStrain(branch);
    for (std::size_t row = 1; row < kappas.size(); ++row)
    {
      SCOPED_TRACE(std::to_string(static_cast<int>(branch)) + ":" + std::to_string(row));
      const auto old = damagedState(table, branch, kappas[row]);
      const double nominal = table.stressValues(branch)[row];
      const double damage = table.response(branch, table.stressAbscissa(branch)[row]).damage.value;
      const double sign = branch == Branch::TENSION ? 1.0 : -1.0;
      const double effective = sign * nominal / (1.0 - damage);
      const auto result = state_integrator.integrate(strainAtEffectiveStress(old, effective), 0.01, old);
      EXPECT_FALSE(result.backbone.plastic);
      EXPECT_NEAR(result.cauchy_stress[0], sign * nominal, nominal * 1.0e-10);
      EXPECT_NEAR(result.cauchy_stress[1], 0.0, nominal * 1.0e-10);
      EXPECT_NEAR(result.cauchy_stress[2], 0.0, nominal * 1.0e-10);
    }
  }
}

TEST(AbaqusCDPEffectiveStrength, QueryPreservesNominalTablesAtKnotsAndBetweenKnots)
{
  const auto table = expertTable();
  for (const auto branch : {Branch::TENSION, Branch::COMPRESSION})
  {
    const auto & kappas = table.equivalentPlasticStrain(branch);
    for (std::size_t row = 0; row < kappas.size(); ++row)
    {
      const auto nominal = table.responseByEquivalentPlasticStrain(branch, kappas[row]);
      const auto effective = table.effectiveStrengthByEquivalentPlasticStrain(branch, kappas[row]);
      EXPECT_DOUBLE_EQ(nominal.stress.value, table.stressValues(branch)[row]);
      EXPECT_NEAR(effective.value * (1.0 - nominal.damage.value), nominal.stress.value,
                  nominal.stress.value * 1.0e-14);
      if (row == 0)
        continue;
      const double middle = 0.5 * (kappas[row - 1] + kappas[row]);
      const auto response = table.responseByEquivalentPlasticStrain(branch, middle);
      const auto q = table.effectiveStrengthByEquivalentPlasticStrain(branch, middle);
      const double expected = response.stress.value / (1.0 - response.damage.value);
      EXPECT_DOUBLE_EQ(q.value, expected);
      EXPECT_DOUBLE_EQ(q.left_derivative, q.right_derivative);
      const double h = (kappas[row] - kappas[row - 1]) * 1.0e-5;
      const double fd =
          (table.effectiveStrengthByEquivalentPlasticStrain(branch, middle + h).value -
           table.effectiveStrengthByEquivalentPlasticStrain(branch, middle - h).value) / (2.0 * h);
      EXPECT_NEAR(q.right_derivative, fd, 1.0e-6 * std::max(1.0, std::abs(fd)));
    }
  }
}

TEST(AbaqusCDPEffectiveStrength, QueryHasCorrectOneSidedKnotAndEndpointDerivatives)
{
  const auto table = expertTable();
  for (const auto branch : {Branch::TENSION, Branch::COMPRESSION})
  {
    const auto & kappas = table.equivalentPlasticStrain(branch);
    for (std::size_t row = 0; row < kappas.size(); ++row)
    {
      SCOPED_TRACE(std::to_string(row));
      const auto q = table.effectiveStrengthByEquivalentPlasticStrain(branch, kappas[row]);
      if (row > 0)
      {
        const double h = (kappas[row] - kappas[row - 1]) * 1.0e-6;
        const double fd = (q.value - table.effectiveStrengthByEquivalentPlasticStrain(
                                        branch, kappas[row] - h).value) / h;
        EXPECT_NEAR(q.left_derivative, fd, 2.0e-5 * std::max(1.0, std::abs(fd)));
      }
      else
        EXPECT_DOUBLE_EQ(q.left_derivative, 0.0);
      if (row + 1 < kappas.size())
      {
        const double h = (kappas[row + 1] - kappas[row]) * 1.0e-6;
        const double fd = (table.effectiveStrengthByEquivalentPlasticStrain(
                               branch, kappas[row] + h).value - q.value) / h;
        EXPECT_NEAR(q.right_derivative, fd, 2.0e-5 * std::max(1.0, std::abs(fd)));
      }
      else
        EXPECT_DOUBLE_EQ(q.right_derivative, 0.0);
    }
    for (const double kappa : {-1.0, 2.0 * kappas.back()})
    {
      const auto q = table.effectiveStrengthByEquivalentPlasticStrain(branch, kappa);
      const auto nominal = table.responseByEquivalentPlasticStrain(branch, kappa);
      EXPECT_DOUBLE_EQ(q.value, nominal.stress.value / (1.0 - nominal.damage.value));
      EXPECT_DOUBLE_EQ(q.left_derivative, 0.0);
      EXPECT_DOUBLE_EQ(q.right_derivative, 0.0);
    }
  }
}

namespace
{
struct ManufacturedReturn
{
  AbaqusCDPStateIntegrator::State old;
  Tensor strain;
  Tensor effective_stress;
  double kappa;
  double nominal_stress;
};

// Construct an exact solution of the backward-Euler plastic equations using
// the public flow potential and sigma/(1-d) as an independent strength oracle.
// The target lies strictly inside a table segment, away from knot derivatives.
ManufacturedReturn
manufacturedReturn(const CDPMaterialTable & table, const Branch branch, const bool smooth = false)
{
  const auto & kappas = table.equivalentPlasticStrain(branch);
  const std::size_t index = kappas.size() / 3;
  const double kappa = 0.5 * (kappas[index] + kappas[index + 1]);
  const auto sample = table.responseByEquivalentPlasticStrain(branch, kappa);
  const double sign = branch == Branch::TENSION ? 1.0 : -1.0;
  const auto inactive_branch = branch == Branch::TENSION ? Branch::COMPRESSION : Branch::TENSION;
  const auto & inactive_kappas = table.equivalentPlasticStrain(inactive_branch);
  const double inactive_kappa = smooth ? 0.5 * (inactive_kappas[3] + inactive_kappas[4]) : 0.0;
  Tensor effective = {sign * sample.stress.value / (1.0 - sample.damage.value),
                      0.0, 0.0, 0.0, 0.0, 0.0};
  if (smooth)
  {
    // A central finite difference is NOT a valid derivative oracle at the two
    // zero principal stresses of exact uniaxial loading (Macaulay/recovery kinks).
    // Use distinct nonzero principals of the same sign to test a smooth branch.
    const auto tension = table.responseByEquivalentPlasticStrain(
        Branch::TENSION, branch == Branch::TENSION ? kappa : inactive_kappa);
    const auto compression = table.responseByEquivalentPlasticStrain(
        Branch::COMPRESSION, branch == Branch::COMPRESSION ? kappa : inactive_kappa);
    const double qt = tension.stress.value / (1.0 - tension.damage.value);
    const double qc = compression.stress.value / (1.0 - compression.damage.value);
    effective = {sign, sign * 0.06, sign * 0.03,
                 sign * 0.01, sign * 0.005, sign * 0.004};
    for (double & value : effective)
      value *= branch == Branch::TENSION ? qt : qc;
    const double homogeneous = AbaqusCDPFormula::yieldFunction(effective, qc, qt, 1.16, 0.667) + qc;
    for (double & value : effective)
      value *= qc / homogeneous;
  }
  const double initial_tension = table.stressValues(Branch::TENSION).front();
  const auto flow = AbaqusCDPFormula::flowPotential(effective, 36.0, 0.1, initial_tension);
  const double multiplier = 0.02 * (kappas[index + 1] - kappas[index]);
  Tensor increment = {};
  for (std::size_t i = 0; i < 6; ++i)
    increment[i] = multiplier * flow.gradient[i] * (i < 3 ? 1.0 : 0.5);
  const auto principal = AbaqusCDPFormula::stressInvariants(increment).principal_stress;
  const double delta_kappa = branch == Branch::TENSION ? principal.back() : -principal.front();
  auto old = damagedState(table, branch, kappa - delta_kappa);
  // Keep the inactive history inside a segment too: at kappa=0 a centered
  // difference averages the constant left extension and the right derivative.
  if (branch == Branch::TENSION)
  {
    old.backbone.compressive_equivalent_plastic_strain = inactive_kappa;
    old.viscous_compression_damage =
        table.responseByEquivalentPlasticStrain(inactive_branch, inactive_kappa).damage.value;
  }
  else
  {
    old.backbone.tensile_equivalent_plastic_strain = inactive_kappa;
    old.viscous_tension_damage =
        table.responseByEquivalentPlasticStrain(inactive_branch, inactive_kappa).damage.value;
  }
  Tensor strain = old.backbone.plastic_strain;
  const double trace = effective[0] + effective[1] + effective[2];
  for (std::size_t i = 0; i < 6; ++i)
    strain[i] += increment[i] + ((1.0 + poisson) * effective[i] -
                                 (i < 3 ? poisson * trace : 0.0)) / modulus;
  return {old, strain, effective, kappa, (1.0 - sample.damage.value) * effective[0]};
}
}

TEST(AbaqusCDPEffectiveStrength, PlasticReturnRecoversManufacturedStateWithAdAndFd)
{
  const auto table = expertTable();
  for (const auto branch : {Branch::TENSION, Branch::COMPRESSION})
    for (const bool automatic : {true, false})
    {
      const auto test = manufacturedReturn(table, branch);
      const AbaqusCDPLocalIntegrator local(table, localParameters(automatic));
      const AbaqusCDPStateIntegrator state_integrator(local, {0.0, 1.0, 0.0, 1.0e-12});
      const auto result = state_integrator.integrate(test.strain, 0.01, test.old);
      EXPECT_TRUE(result.backbone.plastic);
      EXPECT_LT(result.backbone.residual_norm, 1.0e-9);
      EXPECT_NEAR(result.cauchy_stress[0], test.nominal_stress,
                  1.0e-7 * std::abs(test.nominal_stress));
      for (std::size_t i = 0; i < 6; ++i)
        EXPECT_NEAR(result.backbone.effective_stress[i], test.effective_stress[i],
                    1.0e-7 * std::abs(test.effective_stress[0]));
      const double actual_kappa = branch == Branch::TENSION
                                      ? result.state.backbone.tensile_equivalent_plastic_strain
                                      : result.state.backbone.compressive_equivalent_plastic_strain;
      EXPECT_NEAR(actual_kappa, test.kappa, 1.0e-10);
    }
}

TEST(AbaqusCDPEffectiveStrength, DamagedPlasticLocalJacobianMatchesFiniteDifference)
{
  const auto table = expertTable();
  const AbaqusCDPLocalIntegrator local(table, localParameters());
  for (const auto branch : {Branch::TENSION, Branch::COMPRESSION})
  {
    const auto test = manufacturedReturn(table, branch, true);
    const auto matrices = local.localJacobianDiagnostic(test.strain, test.old.backbone);
    for (std::size_t row = 0; row < local.local_size; ++row)
      for (std::size_t col = 0; col < local.local_size; ++col)
      {
        SCOPED_TRACE(std::to_string(static_cast<int>(branch)) + ":" +
                     std::to_string(row) + ":" + std::to_string(col));
        EXPECT_NEAR(matrices.automatic_differentiation[row][col],
                    matrices.finite_difference[row][col],
                    5.0e-4 * std::max(1.0, std::abs(matrices.finite_difference[row][col])));
      }
  }
}

TEST(AbaqusCDPEffectiveStrength, DamagedHistoryAndStrainTransitionMatchFiniteDifference)
{
  const auto table = expertTable();
  const AbaqusCDPLocalIntegrator local(table, localParameters());
  for (const auto branch : {Branch::TENSION, Branch::COMPRESSION})
  {
    const auto test = manufacturedReturn(table, branch, true);
    const auto linearized = local.integrateLinearized(test.strain, test.old.backbone);
    ASSERT_TRUE(linearized.result.plastic);
    // Include both old kappas, not just six strain columns: these direct
    // contributions are needed when chaining the uniform substep tangent.
    for (std::size_t col = 0; col < local.transition_size; ++col)
    {
      SCOPED_TRACE(std::to_string(static_cast<int>(branch)) + ":" + std::to_string(col));
      auto plus_strain = test.strain;
      auto minus_strain = test.strain;
      auto plus_old = test.old.backbone;
      auto minus_old = test.old.backbone;
      const double h = 1.0e-10;
      if (col < 6)
      {
        plus_strain[col] += h;
        minus_strain[col] -= h;
      }
      else if (col < 12)
      {
        plus_old.plastic_strain[col - 6] += h;
        minus_old.plastic_strain[col - 6] -= h;
      }
      else if (col == 12)
      {
        if (minus_old.tensile_equivalent_plastic_strain == 0.0)
          continue; // inactive history at its nonnegative boundary
        plus_old.tensile_equivalent_plastic_strain += h;
        minus_old.tensile_equivalent_plastic_strain -= h;
      }
      else
      {
        if (minus_old.compressive_equivalent_plastic_strain == 0.0)
          continue;
        plus_old.compressive_equivalent_plastic_strain += h;
        minus_old.compressive_equivalent_plastic_strain -= h;
      }
      const auto plus = local.integrate(plus_strain, plus_old);
      const auto minus = local.integrate(minus_strain, minus_old);
      for (std::size_t row = 0; row < 6; ++row)
      {
        const double fd = (plus.effective_stress[row] - minus.effective_stress[row]) / (2.0 * h);
        EXPECT_NEAR(linearized.derivative[col][row], fd, 2.0e-5 * modulus);
        const double plastic_fd = (plus.state.plastic_strain[row] -
                                    minus.state.plastic_strain[row]) / (2.0 * h);
        EXPECT_NEAR(linearized.derivative[col][6 + row], plastic_fd, 2.0e-5);
      }
      EXPECT_NEAR(linearized.derivative[col][12],
                  (plus.state.tensile_equivalent_plastic_strain -
                   minus.state.tensile_equivalent_plastic_strain) / (2.0 * h), 2.0e-5);
      EXPECT_NEAR(linearized.derivative[col][13],
                  (plus.state.compressive_equivalent_plastic_strain -
                   minus.state.compressive_equivalent_plastic_strain) / (2.0 * h), 2.0e-5);
    }
  }
}

TEST(AbaqusCDPEffectiveStrength, DamagedCauchyTangentMatchesFiniteDifferenceWithAndWithoutViscosity)
{
  const auto table = expertTable();
  const AbaqusCDPLocalIntegrator local(table, localParameters());
  for (const auto branch : {Branch::TENSION, Branch::COMPRESSION})
    for (const double viscosity : {0.0, 0.0005})
    {
      const auto test = manufacturedReturn(table, branch, true);
      const AbaqusCDPStateIntegrator integrator(local, {0.0, 1.0, viscosity, 1.0e-12});
      const auto linearized = integrator.integrateLinearized(test.strain, 0.01, test.old);
      ASSERT_TRUE(linearized.result.backbone.plastic);
      for (std::size_t col = 0; col < 6; ++col)
      {
        auto plus = test.strain;
        auto minus = test.strain;
        const double h = 1.0e-10;
        plus[col] += h;
        minus[col] -= h;
        const auto plus_result = integrator.integrate(plus, 0.01, test.old);
        const auto minus_result = integrator.integrate(minus, 0.01, test.old);
        for (std::size_t row = 0; row < 6; ++row)
        {
          const double fd = (plus_result.cauchy_stress[row] - minus_result.cauchy_stress[row]) /
                              (2.0 * h);
          EXPECT_NEAR(linearized.derivative[col][row], fd, 2.0e-5 * modulus);
        }
      }
    }
}
