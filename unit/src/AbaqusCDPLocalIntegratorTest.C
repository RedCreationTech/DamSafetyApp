#include "AbaqusCDPLocalIntegrator.h"

#include "gtest/gtest.h"

#include <cmath>
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
  const std::array<AbaqusCDPLocalIntegrator::SymmetricTensor, 3> strains = {
      uniaxialElasticStrain(1.05 * initial_tension),
      uniaxialElasticStrain(-1.05 * initial_compression),
      AbaqusCDPLocalIntegrator::SymmetricTensor{2.0e-4,
                                               -8.0e-5,
                                               -4.0e-5,
                                               3.0e-5,
                                               -1.0e-5,
                                               2.0e-5}};

  for (const auto & strain : strains)
  {
    const auto automatic_result = automatic.integrateLinearized(strain, {});
    const auto reference_result = finite_difference.integrateLinearized(strain, {});
    ASSERT_EQ(automatic_result.result.active_branch, reference_result.result.active_branch);
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
    for (std::size_t input = 0; input < AbaqusCDPLocalIntegrator::transition_size; ++input)
      for (std::size_t output = 0; output < AbaqusCDPLocalIntegrator::transition_size; ++output)
        EXPECT_NEAR(automatic_result.derivative[input][output],
                    reference_result.derivative[input][output],
                    5.0e-3 *
                        std::max(1.0, std::abs(reference_result.derivative[input][output])));
  }
}
