#include "AbaqusCDPStateIntegrator.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace
{
constexpr double youngs_modulus = 3.04e10;
constexpr double poissons_ratio = 0.2;
const std::string data_dir = "test/tests/cdp_material_table/data/";

CDPMaterialTable
stateReferenceTable()
{
  return CDPMaterialTable(data_dir + "compression_hardening.csv",
                          data_dir + "compression_damage.csv",
                          data_dir + "tension_stiffening.csv",
                          data_dir + "tension_damage.csv",
                          youngs_modulus);
}

AbaqusCDPLocalIntegrator::Parameters
localParameters(const unsigned int maximum_iterations = 40)
{
  return {youngs_modulus,
          poissons_ratio,
          36.31,
          0.1,
          1.16,
          0.667,
          maximum_iterations,
          1.0e-9,
          1.0e-7,
          1.0e-6};
}

AbaqusCDPStateIntegrator::Parameters
stateParameters(const double relaxation_time,
                const double tension_recovery = 1.0,
                const double compression_recovery = 0.0)
{
  return {tension_recovery, compression_recovery, relaxation_time, 1.0e-12};
}

AbaqusCDPStateIntegrator::SymmetricTensor
stateUniaxialElasticStrain(const double stress)
{
  return {stress / youngs_modulus,
          -poissons_ratio * stress / youngs_modulus,
          -poissons_ratio * stress / youngs_modulus,
          0.0,
          0.0,
          0.0};
}

AbaqusCDPStateIntegrator::SymmetricTensor
addStateStrain(const AbaqusCDPStateIntegrator::SymmetricTensor & left,
               const AbaqusCDPStateIntegrator::SymmetricTensor & right)
{
  AbaqusCDPStateIntegrator::SymmetricTensor result;
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = left[i] + right[i];
  return result;
}
}

TEST(AbaqusCDPStateIntegrator, ZeroViscosityMatchesBackboneAndAppliesDamage)
{
  const auto table = stateReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, localParameters());
  const AbaqusCDPStateIntegrator integrator(local, stateParameters(0.0));
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;

  const auto result =
      integrator.integrate(stateUniaxialElasticStrain(1.05 * initial_tension), 1.0e-3, {});

  EXPECT_TRUE(result.backbone.plastic);
  EXPECT_EQ(result.state.viscous_plastic_strain, result.backbone.state.plastic_strain);
  EXPECT_DOUBLE_EQ(result.state.viscous_tension_damage,
                   result.backbone.backbone_tension_damage);
  EXPECT_DOUBLE_EQ(result.state.viscous_compression_damage,
                   result.backbone.backbone_compression_damage);
  for (std::size_t i = 0; i < result.viscous_effective_stress.size(); ++i)
    EXPECT_NEAR(result.viscous_effective_stress[i],
                result.backbone.effective_stress[i],
                1.0e-12 * initial_tension);
  EXPECT_GT(result.state.viscous_tension_damage, 0.0);
  for (std::size_t i = 0; i < result.cauchy_stress.size(); ++i)
    EXPECT_NEAR(result.cauchy_stress[i],
                result.damage.stiffness_factor * result.viscous_effective_stress[i],
                1.0e-12 * std::max(1.0, std::abs(result.viscous_effective_stress[i])));
  EXPECT_TRUE(std::isinf(result.dt_over_relaxation_time));
  EXPECT_DOUBLE_EQ(result.plastic_strain_lag_norm, 0.0);
}

TEST(AbaqusCDPStateIntegrator, AppliesBackwardEulerToPlasticStrainAndBothDamageBranches)
{
  const auto table = stateReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, localParameters());
  const double mu = 5.0e-4;
  const double dt = 5.0e-4;
  const AbaqusCDPStateIntegrator integrator(local, stateParameters(mu));
  const double initial_compression =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::COMPRESSION, 0.0)
          .stress.value;

  const auto result =
      integrator.integrate(stateUniaxialElasticStrain(-1.05 * initial_compression), dt, {});

  for (std::size_t i = 0; i < result.state.viscous_plastic_strain.size(); ++i)
    EXPECT_NEAR(result.state.viscous_plastic_strain[i],
                0.5 * result.backbone.state.plastic_strain[i],
                1.0e-14);
  EXPECT_NEAR(result.state.viscous_compression_damage,
              0.5 * result.backbone.backbone_compression_damage,
              1.0e-14);
  EXPECT_NEAR(result.state.viscous_tension_damage,
              0.5 * result.backbone.backbone_tension_damage,
              1.0e-14);
  EXPECT_DOUBLE_EQ(result.dt_over_relaxation_time, 1.0);
  EXPECT_GT(result.plastic_strain_lag_norm, 0.0);
  EXPECT_GE(result.compression_damage_lag, 0.0);
}

TEST(AbaqusCDPStateIntegrator, ZeroTimeStepKeepsViscousHistory)
{
  const auto table = stateReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, localParameters());
  const AbaqusCDPStateIntegrator integrator(local, stateParameters(5.0e-4));
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto first =
      integrator.integrate(stateUniaxialElasticStrain(1.05 * initial_tension), 2.5e-4, {});
  const auto hold = integrator.integrate(
      stateUniaxialElasticStrain(1.05 * initial_tension), 0.0, first.state);

  EXPECT_EQ(hold.state.viscous_plastic_strain, first.state.viscous_plastic_strain);
  EXPECT_DOUBLE_EQ(hold.state.viscous_tension_damage, first.state.viscous_tension_damage);
  EXPECT_DOUBLE_EQ(hold.state.viscous_compression_damage,
                   first.state.viscous_compression_damage);
  EXPECT_DOUBLE_EQ(hold.dt_over_relaxation_time, 0.0);
}

TEST(AbaqusCDPStateIntegrator, TensionRecoveryControlsPreviouslyAccumulatedCompressionDamage)
{
  const auto table = stateReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, localParameters());
  const AbaqusCDPStateIntegrator accumulating(local, stateParameters(0.0));
  const double initial_compression =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::COMPRESSION, 0.0)
          .stress.value;
  const auto compressed = accumulating.integrate(
      stateUniaxialElasticStrain(-1.05 * initial_compression), 1.0e-3, {});
  ASSERT_GT(compressed.state.viscous_compression_damage, 0.0);

  const auto reverse_tension = addStateStrain(
      compressed.state.backbone.plastic_strain, stateUniaxialElasticStrain(1.0e6));
  const AbaqusCDPStateIntegrator full_recovery(local, stateParameters(0.0, 1.0, 0.0));
  const AbaqusCDPStateIntegrator no_recovery(local, stateParameters(0.0, 0.0, 0.0));
  const auto recovered = full_recovery.integrate(reverse_tension, 1.0e-3, compressed.state);
  const auto unrecovered = no_recovery.integrate(reverse_tension, 1.0e-3, compressed.state);

  EXPECT_NEAR(recovered.damage.tension_weight, 1.0, 1.0e-12);
  EXPECT_NEAR(recovered.damage.stiffness_factor, 1.0, 1.0e-12);
  EXPECT_NEAR(unrecovered.damage.stiffness_factor,
              1.0 - compressed.state.viscous_compression_damage,
              1.0e-12);
  EXPECT_GT(recovered.cauchy_stress[0], unrecovered.cauchy_stress[0]);
}

TEST(AbaqusCDPStateIntegrator, FailedBackboneSolveDoesNotMutateOldState)
{
  const auto table = stateReferenceTable();
  const AbaqusCDPLocalIntegrator failing_local(table, localParameters(0));
  const AbaqusCDPStateIntegrator integrator(failing_local, stateParameters(5.0e-4));
  AbaqusCDPStateIntegrator::State old_state;
  const auto snapshot = old_state;
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;

  EXPECT_THROW(integrator.integrate(
                   stateUniaxialElasticStrain(1.05 * initial_tension), 1.0e-4, old_state),
               std::runtime_error);
  EXPECT_EQ(old_state.backbone.plastic_strain, snapshot.backbone.plastic_strain);
  EXPECT_EQ(old_state.viscous_plastic_strain, snapshot.viscous_plastic_strain);
  EXPECT_DOUBLE_EQ(old_state.viscous_tension_damage, snapshot.viscous_tension_damage);
}

TEST(AbaqusCDPStateIntegrator, RejectsInvalidStateAndParameters)
{
  const auto table = stateReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, localParameters());
  EXPECT_THROW(AbaqusCDPStateIntegrator(local, stateParameters(-1.0)), std::runtime_error);

  const AbaqusCDPStateIntegrator integrator(local, stateParameters(5.0e-4));
  AbaqusCDPStateIntegrator::State invalid;
  invalid.viscous_compression_damage = 1.0;
  EXPECT_THROW(integrator.integrate({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 1.0e-4, invalid),
               std::runtime_error);
  EXPECT_THROW(integrator.integrate({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, -1.0, {}),
               std::runtime_error);
}
