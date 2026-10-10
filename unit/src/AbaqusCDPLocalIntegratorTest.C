#include "AbaqusCDPLocalIntegrator.h"

#include "gtest/gtest.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>

namespace
{
constexpr double youngs_modulus = 3.04e10;
constexpr double poissons_ratio = 0.2;
const std::string data_dir = "test/tests/cdp_material_table/data/";

CDPMaterialTable
referenceTable()
{
  return CDPMaterialTable(data_dir + "compression_hardening.csv",
                          data_dir + "compression_damage.csv",
                          data_dir + "tension_stiffening.csv",
                          data_dir + "tension_damage.csv",
                          youngs_modulus);
}

AbaqusCDPLocalIntegrator::Parameters
parameters(const unsigned int maximum_iterations = 40,
           const bool use_automatic_differentiation_jacobian = true)
{
  auto result = AbaqusCDPLocalIntegrator::Parameters{youngs_modulus,
                                                     poissons_ratio,
                                                     36.31,
                                                     0.1,
                                                     1.16,
                                                     0.667,
                                                     maximum_iterations,
                                                     1.0e-9,
                                                     1.0e-7,
                                                     1.0e-6};
  result.use_automatic_differentiation_jacobian = use_automatic_differentiation_jacobian;
  return result;
}

AbaqusCDPLocalIntegrator::SymmetricTensor
uniaxialElasticStrain(const double stress)
{
  return {stress / youngs_modulus,
          -poissons_ratio * stress / youngs_modulus,
          -poissons_ratio * stress / youngs_modulus,
          0.0,
          0.0,
          0.0};
}
}

TEST(AbaqusCDPLocalIntegrator, ReturnsExactElasticTrialAndUnchangedState)
{
  const auto table = referenceTable();
  const AbaqusCDPLocalIntegrator integrator(table, parameters());
  AbaqusCDPLocalIntegrator::State old_state;
  const auto result = integrator.integrate(uniaxialElasticStrain(1.0e6), old_state);

  EXPECT_FALSE(result.plastic);
  EXPECT_EQ(result.active_branch, AbaqusCDPLocalIntegrator::ActiveBranch::ELASTIC);
  EXPECT_DOUBLE_EQ(result.effective_stress[0], 1.0e6);
  EXPECT_NEAR(result.effective_stress[1], 0.0, 1.0e-10);
  EXPECT_NEAR(result.effective_stress[2], 0.0, 1.0e-10);
  EXPECT_EQ(result.state.plastic_strain, old_state.plastic_strain);
  EXPECT_DOUBLE_EQ(result.state.tensile_equivalent_plastic_strain, 0.0);
  EXPECT_DOUBLE_EQ(result.state.compressive_equivalent_plastic_strain, 0.0);
}

TEST(AbaqusCDPLocalIntegrator, ConvergesTensilePlasticCorrectionAndHardeningResidual)
{
  const auto table = referenceTable();
  const AbaqusCDPLocalIntegrator integrator(table, parameters());
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto result =
      integrator.integrate(uniaxialElasticStrain(1.05 * initial_tension), {});

  EXPECT_TRUE(result.plastic);
  EXPECT_EQ(result.active_branch, AbaqusCDPLocalIntegrator::ActiveBranch::TENSION);
  EXPECT_GT(result.plastic_multiplier, 0.0);
  EXPECT_GT(result.state.tensile_equivalent_plastic_strain, 0.0);
  EXPECT_NEAR(result.state.compressive_equivalent_plastic_strain, 0.0, 1.0e-14);
  EXPECT_LT(result.residual_norm, 1.0e-9);
  EXPECT_NEAR(result.final_yield / initial_tension, 0.0, 1.0e-8);
}

TEST(AbaqusCDPLocalIntegrator, ConvergesCompressivePlasticCorrectionAndHardeningResidual)
{
  const auto table = referenceTable();
  const AbaqusCDPLocalIntegrator integrator(table, parameters());
  const double initial_compression = table
                                         .responseByEquivalentPlasticStrain(
                                             CDPMaterialTable::Branch::COMPRESSION, 0.0)
                                         .stress.value;
  const auto result =
      integrator.integrate(uniaxialElasticStrain(-1.05 * initial_compression), {});

  EXPECT_TRUE(result.plastic);
  EXPECT_EQ(result.active_branch, AbaqusCDPLocalIntegrator::ActiveBranch::COMPRESSION);
  EXPECT_GT(result.plastic_multiplier, 0.0);
  EXPECT_GT(result.state.compressive_equivalent_plastic_strain, 0.0);
  EXPECT_NEAR(result.state.tensile_equivalent_plastic_strain, 0.0, 1.0e-14);
  EXPECT_LT(result.residual_norm, 1.0e-9);
  EXPECT_NEAR(result.final_yield / initial_compression, 0.0, 1.0e-8);
}

TEST(AbaqusCDPLocalIntegrator, FailedSolveDoesNotMutateOldStateAndRetryIsDeterministic)
{
  const auto table = referenceTable();
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto strain = uniaxialElasticStrain(1.05 * initial_tension);
  AbaqusCDPLocalIntegrator::State old_state;
  const auto snapshot = old_state;

  const AbaqusCDPLocalIntegrator failing(table, parameters(0));
  try
  {
    failing.integrate(strain, old_state);
    FAIL() << "expected local integration failure";
  }
  catch (const std::runtime_error & error)
  {
    const std::string diagnostic = error.what();
    EXPECT_NE(diagnostic.find("residual="), std::string::npos);
    EXPECT_NE(diagnostic.find("plastic_multiplier="), std::string::npos);
    EXPECT_NE(diagnostic.find("kappa_t="), std::string::npos);
    EXPECT_NE(diagnostic.find("kappa_c="), std::string::npos);
    EXPECT_NE(diagnostic.find("branch="), std::string::npos);
  }
  EXPECT_EQ(old_state.plastic_strain, snapshot.plastic_strain);
  EXPECT_DOUBLE_EQ(old_state.tensile_equivalent_plastic_strain,
                   snapshot.tensile_equivalent_plastic_strain);
  EXPECT_DOUBLE_EQ(old_state.compressive_equivalent_plastic_strain,
                   snapshot.compressive_equivalent_plastic_strain);

  const AbaqusCDPLocalIntegrator working(table, parameters());
  const auto first = working.integrate(strain, old_state);
  const auto second = working.integrate(strain, old_state);
  EXPECT_EQ(first.effective_stress, second.effective_stress);
  EXPECT_EQ(first.state.plastic_strain, second.state.plastic_strain);
  EXPECT_DOUBLE_EQ(first.state.tensile_equivalent_plastic_strain,
                   second.state.tensile_equivalent_plastic_strain);
  EXPECT_DOUBLE_EQ(first.state.compressive_equivalent_plastic_strain,
                   second.state.compressive_equivalent_plastic_strain);
}

TEST(AbaqusCDPLocalIntegrator, RejectsInvalidOldStateBeforeTrialEvaluation)
{
  const auto table = referenceTable();
  const AbaqusCDPLocalIntegrator integrator(table, parameters());
  AbaqusCDPLocalIntegrator::State invalid;
  invalid.tensile_equivalent_plastic_strain = -1.0;
  EXPECT_THROW(integrator.integrate({0, 0, 0, 0, 0, 0}, invalid), std::runtime_error);
}

TEST(AbaqusCDPLocalIntegrator, AutomaticDifferentiationJacobianMatchesFiniteDifference)
{
  const auto table = referenceTable();
  const AbaqusCDPLocalIntegrator integrator(table, parameters());
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  auto strain = uniaxialElasticStrain(1.08 * initial_tension);
  strain[3] = 0.08 * initial_tension / youngs_modulus;

  const auto diagnostic = integrator.localJacobianDiagnostic(strain, {});
  double maximum_absolute_error = 0.0;
  double maximum_scaled_error = 0.0;
  for (std::size_t row = 0; row < AbaqusCDPLocalIntegrator::local_size; ++row)
    for (std::size_t column = 0; column < AbaqusCDPLocalIntegrator::local_size; ++column)
    {
      const double automatic = diagnostic.automatic_differentiation[row][column];
      const double reference = diagnostic.finite_difference[row][column];
      const double difference = std::abs(automatic - reference);
      maximum_absolute_error = std::max(maximum_absolute_error, difference);
      maximum_scaled_error =
          std::max(maximum_scaled_error, difference / std::max(1.0, std::abs(reference)));
    }
  EXPECT_LE(maximum_absolute_error, 1.0e-6);
  EXPECT_LE(maximum_scaled_error, 5.0e-4);
}

TEST(AbaqusCDPLocalIntegrator, AutomaticDifferentiationAndFiniteDifferenceReachSameRoots)
{
  const auto table = referenceTable();
  const AbaqusCDPLocalIntegrator automatic(table, parameters(40, true));
  const AbaqusCDPLocalIntegrator finite_difference(table, parameters(40, false));
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const double initial_compression = table
                                         .responseByEquivalentPlasticStrain(
                                             CDPMaterialTable::Branch::COMPRESSION, 0.0)
                                         .stress.value;
  auto sheared_tension = uniaxialElasticStrain(1.08 * initial_tension);
  sheared_tension[3] = 0.08 * initial_tension / youngs_modulus;
  const std::array<AbaqusCDPLocalIntegrator::SymmetricTensor, 3> strains = {
      uniaxialElasticStrain(1.05 * initial_tension),
      uniaxialElasticStrain(-1.05 * initial_compression),
      sheared_tension};

  for (const auto & strain : strains)
  {
    const auto automatic_result = automatic.integrateLinearized(strain, {});
    const auto reference_result = finite_difference.integrateLinearized(strain, {});
    ASSERT_EQ(automatic_result.result.active_branch, reference_result.result.active_branch);
    EXPECT_GT(automatic_result.result.local_factorizations, 0u);
    EXPECT_EQ(automatic_result.result.local_backsolves,
              automatic_result.result.local_factorizations + 7u);
    EXPECT_GT(reference_result.result.local_factorizations, 0u);
    EXPECT_EQ(reference_result.result.local_backsolves,
              reference_result.result.local_factorizations + 7u);
    for (std::size_t i = 0; i < 6; ++i)
    {
      EXPECT_NEAR(automatic_result.result.effective_stress[i],
                  reference_result.result.effective_stress[i],
                  1.0e-7 * std::max(1.0, std::abs(reference_result.result.effective_stress[i])));
      EXPECT_NEAR(automatic_result.result.state.plastic_strain[i],
                  reference_result.result.state.plastic_strain[i],
                  1.0e-10);
    }
    EXPECT_NEAR(automatic_result.result.state.tensile_equivalent_plastic_strain,
                reference_result.result.state.tensile_equivalent_plastic_strain,
                1.0e-10);
    EXPECT_NEAR(automatic_result.result.state.compressive_equivalent_plastic_strain,
                reference_result.result.state.compressive_equivalent_plastic_strain,
                1.0e-10);
  }
}

TEST(AbaqusCDPLocalIntegrator, ReusedPlasticStrainColumnsPreserveTransitionAndEightBacksolves)
{
  const auto table = referenceTable();
  for (const bool automatic : {true, false})
  {
    const AbaqusCDPLocalIntegrator integrator(table, parameters(40, automatic));
    const double tension = table.responseByEquivalentPlasticStrain(
        CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
    const double compression = table.responseByEquivalentPlasticStrain(
        CDPMaterialTable::Branch::COMPRESSION, 0.0).stress.value;
    auto mixed = uniaxialElasticStrain(1.08 * tension);
    mixed[3] = 0.08 * tension / youngs_modulus;
    auto smooth_compression = uniaxialElasticStrain(-1.05 * compression);
    // Split the repeated transverse principal stresses: centered differences at
    // the exact compression meridian need not equal a selected branch derivative.
    smooth_compression[3] = 0.07 * compression / youngs_modulus;
    smooth_compression[4] = -0.03 * compression / youngs_modulus;
    for (const auto & strain : {uniaxialElasticStrain(1.05 * tension),
                               smooth_compression, mixed})
    {
      const auto plain = integrator.integrate(strain, {});
      const auto result = integrator.integrateLinearized(strain, {});
      ASSERT_TRUE(result.result.plastic);
      EXPECT_EQ(result.result.local_backsolves - plain.local_backsolves, 8u);
      EXPECT_EQ(result.result.effective_stress, plain.effective_stress);
      EXPECT_EQ(result.result.state.plastic_strain, plain.state.plastic_strain);
      for (std::size_t column = 0; column < 6; ++column)
      {
        auto plus = AbaqusCDPLocalIntegrator::State{};
        auto minus = plus;
        constexpr double h = 1.0e-9;
        plus.plastic_strain[column] += h;
        minus.plastic_strain[column] -= h;
        const auto rp = integrator.integrate(strain, plus);
        const auto rm = integrator.integrate(strain, minus);
        ASSERT_EQ(rp.active_branch, result.result.active_branch);
        ASSERT_EQ(rm.active_branch, result.result.active_branch);
        for (std::size_t row = 0; row < 6; ++row)
        {
          EXPECT_DOUBLE_EQ(result.derivative[6 + column][row],
                           -result.derivative[column][row]);
          EXPECT_NEAR(result.derivative[6 + column][6 + row] +
                          result.derivative[column][6 + row],
                      column == row ? 1.0 : 0.0, 1.0e-14);
          EXPECT_NEAR(result.derivative[6 + column][row],
                      (rp.effective_stress[row] - rm.effective_stress[row]) / (2 * h),
                      2.0e-4 * youngs_modulus);
          EXPECT_NEAR(result.derivative[6 + column][6 + row],
                      (rp.state.plastic_strain[row] - rm.state.plastic_strain[row]) / (2 * h),
                      2.0e-4);
        }
        EXPECT_NEAR(result.derivative[6 + column][12],
                    (rp.state.tensile_equivalent_plastic_strain -
                     rm.state.tensile_equivalent_plastic_strain) / (2 * h), 2.0e-4);
        EXPECT_NEAR(result.derivative[6 + column][13],
                    (rp.state.compressive_equivalent_plastic_strain -
                     rm.state.compressive_equivalent_plastic_strain) / (2 * h), 2.0e-4);
      }
    }
  }
}

TEST(AbaqusCDPLocalIntegrator, PassiveNewtonTracePreservesRootAndFullTransition)
{
  const auto table = referenceTable();
  for (const bool automatic : {true, false})
  {
    const AbaqusCDPLocalIntegrator integrator(table, parameters(40, automatic));
    for (const auto branch : {CDPMaterialTable::Branch::TENSION,
                              CDPMaterialTable::Branch::COMPRESSION})
    {
      const double strength = table.responseByEquivalentPlasticStrain(branch, 0.0).stress.value;
      const auto target = uniaxialElasticStrain(
          (branch == CDPMaterialTable::Branch::TENSION ? 1.05 : -1.05) * strength);
      AbaqusCDPLocalIntegrator::State old;
      const auto snapshot = old;
      const auto control = integrator.integrateLinearized(target, old);
      AbaqusCDPLocalIntegrator::NewtonTrace trace;
      const auto observed = integrator.integrateLinearized(target, old, &trace);
      EXPECT_EQ(control.derivative, observed.derivative);
      EXPECT_EQ(control.result.effective_stress, observed.result.effective_stress);
      EXPECT_EQ(control.result.state.plastic_strain, observed.result.state.plastic_strain);
      EXPECT_EQ(control.result.state.tensile_equivalent_plastic_strain,
                observed.result.state.tensile_equivalent_plastic_strain);
      EXPECT_EQ(control.result.state.compressive_equivalent_plastic_strain,
                observed.result.state.compressive_equivalent_plastic_strain);
      EXPECT_EQ(control.result.residual_norm, observed.result.residual_norm);
      EXPECT_EQ(control.result.iterations, observed.result.iterations);
      EXPECT_EQ(control.result.jacobian_fallbacks, observed.result.jacobian_fallbacks);
      EXPECT_EQ(control.result.local_factorizations, observed.result.local_factorizations);
      EXPECT_EQ(control.result.local_backsolves, observed.result.local_backsolves);
      EXPECT_EQ(control.result.automatic_jacobian_evaluations,
                observed.result.automatic_jacobian_evaluations);
      EXPECT_EQ(control.result.finite_difference_jacobian_evaluations,
                observed.result.finite_difference_jacobian_evaluations);
      EXPECT_EQ(old.plastic_strain, snapshot.plastic_strain);
      EXPECT_EQ(old.tensile_equivalent_plastic_strain, snapshot.tensile_equivalent_plastic_strain);
      EXPECT_EQ(old.compressive_equivalent_plastic_strain, snapshot.compressive_equivalent_plastic_strain);
      ASSERT_FALSE(trace.empty());
      EXPECT_EQ(trace.front().event, "initial");
      EXPECT_EQ(trace.back().event, "tangent_jacobian");
      bool direction_seen = false;
      for (const auto & entry : trace)
      {
        if (entry.event != "direction") continue;
        direction_seen = true;
        ASSERT_TRUE(entry.has_jacobian);
        ASSERT_TRUE(entry.has_direction);
        // The recorded direction must solve the actual recorded nine-equation system.
        for (std::size_t row = 0; row < 9; ++row)
        {
          double value = entry.residual[row], magnitude = std::abs(value);
          for (std::size_t column = 0; column < 9; ++column)
          {
            const double term = entry.jacobian[row][column] * entry.direction[column];
            value += term;
            magnitude += std::abs(term);
          }
          EXPECT_LE(std::abs(value), 1e-12 * std::max(1.0, magnitude));
        }
      }
      EXPECT_TRUE(direction_seen);
    }
  }
}

TEST(AbaqusCDPLocalIntegrator, PassiveNewtonTraceKeepsFailureAndImmutableHistory)
{
  const auto table = referenceTable();
  const AbaqusCDPLocalIntegrator integrator(table, parameters(0));
  const double strength = table.responseByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto target = uniaxialElasticStrain(1.05 * strength);
  AbaqusCDPLocalIntegrator::State old;
  const auto snapshot = old;
  std::string control_error, observed_error;
  try { integrator.integrateLinearized(target, old); }
  catch (const std::runtime_error & error) { control_error = error.what(); }
  AbaqusCDPLocalIntegrator::NewtonTrace trace;
  try { integrator.integrateLinearized(target, old, &trace); }
  catch (const std::runtime_error & error) { observed_error = error.what(); }
  ASSERT_FALSE(control_error.empty());
  EXPECT_EQ(observed_error, control_error);
  ASSERT_EQ(trace.size(), 2u);
  EXPECT_EQ(trace.front().event, "initial");
  EXPECT_EQ(trace.back().event, "iteration_limit");
  EXPECT_NE(control_error.find(trace.back().error), std::string::npos);
  EXPECT_TRUE(trace.back().evaluated);
  EXPECT_FALSE(trace.back().accepted);
  EXPECT_EQ(old.plastic_strain, snapshot.plastic_strain);
  EXPECT_EQ(old.tensile_equivalent_plastic_strain, snapshot.tensile_equivalent_plastic_strain);
  EXPECT_EQ(old.compressive_equivalent_plastic_strain, snapshot.compressive_equivalent_plastic_strain);
}
