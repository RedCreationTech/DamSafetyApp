#include "CDPDiagnostics.h"
#include "AbaqusCDPSubstepIntegrator.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>

namespace
{
constexpr double substep_youngs_modulus = 3.04e10;
constexpr double substep_poissons_ratio = 0.2;
const std::string substep_data_dir = "test/tests/cdp_material_table/data/";
const std::string captured_history_data_dir =
    "test/tests/cdp_material_table/uniaxial_tension_20260830/";

CDPMaterialTable
substepReferenceTable()
{
  return CDPMaterialTable(substep_data_dir + "compression_hardening.csv",
                          substep_data_dir + "compression_damage.csv",
                          substep_data_dir + "tension_stiffening.csv",
                          substep_data_dir + "tension_damage.csv",
                          substep_youngs_modulus);
}

CDPMaterialTable
capturedHistoryTable()
{
  return CDPMaterialTable(captured_history_data_dir + "compression_hardening.csv",
                          captured_history_data_dir + "compression_damage.csv",
                          captured_history_data_dir + "tension_stiffening.csv",
                          captured_history_data_dir + "tension_damage.csv",
                          2.97915e10);
}

AbaqusCDPLocalIntegrator::Parameters
substepLocalParameters(const unsigned int maximum_iterations = 40)
{
  return {substep_youngs_modulus,
          substep_poissons_ratio,
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
substepStateParameters(const double relaxation_time = 0.0)
{
  return {1.0, 0.0, relaxation_time, 1.0e-12};
}

AbaqusCDPSubstepIntegrator::SymmetricTensor
substepUniaxialElasticStrain(const double stress)
{
  return {stress / substep_youngs_modulus,
          -substep_poissons_ratio * stress / substep_youngs_modulus,
          -substep_poissons_ratio * stress / substep_youngs_modulus,
          0.0,
          0.0,
          0.0};
}

AbaqusCDPSubstepIntegrator::SymmetricTensor
substepScale(const AbaqusCDPSubstepIntegrator::SymmetricTensor & tensor, const double factor)
{
  AbaqusCDPSubstepIntegrator::SymmetricTensor result;
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = factor * tensor[i];
  return result;
}

AbaqusCDPSubstepIntegrator::SymmetricTensor
substepDifference(const AbaqusCDPSubstepIntegrator::SymmetricTensor & left,
                  const AbaqusCDPSubstepIntegrator::SymmetricTensor & right)
{
  AbaqusCDPSubstepIntegrator::SymmetricTensor result;
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = left[i] - right[i];
  return result;
}
}

TEST(AbaqusCDPSubstepIntegrator, ProactiveBinarySubstepsMatchManualSequentialIntegration)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, substepStateParameters(5.0e-4));
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto target = substepUniaxialElasticStrain(1.05 * initial_tension);
  const AbaqusCDPSubstepIntegrator integrator(
      state_integrator, {16, std::abs(target[0]) / 3.0, 1.0e-8});

  const auto result = integrator.integrate({}, target, 1.0e-3, {});
  EXPECT_EQ(result.accepted_substeps, 4u);
  EXPECT_EQ(result.cutback_count, 0u);
  EXPECT_EQ(result.attempted_partitions, 1u);
  EXPECT_TRUE(result.proactively_partitioned);

  AbaqusCDPStateIntegrator::State manual_state;
  std::optional<AbaqusCDPStateIntegrator::Result> manual_result;
  for (unsigned int i = 1; i <= 4; ++i)
  {
    manual_result = state_integrator.integrate(
        substepScale(target, static_cast<double>(i) / 4.0), 2.5e-4, manual_state);
    manual_state = manual_result->state;
  }
  ASSERT_TRUE(manual_result.has_value());
  EXPECT_EQ(result.final_result.cauchy_stress, manual_result->cauchy_stress);
  EXPECT_EQ(result.final_result.state.viscous_plastic_strain,
            manual_result->state.viscous_plastic_strain);
  EXPECT_DOUBLE_EQ(result.final_result.state.viscous_tension_damage,
                   manual_result->state.viscous_tension_damage);
}

TEST(AbaqusCDPSubstepIntegrator, DeferredViscousUpdateCommitsOnceAtGlobalStep)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, substepStateParameters(5.0e-4));
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto target = substepUniaxialElasticStrain(1.05 * initial_tension);
  const double strain_limit = std::abs(target[0]) / 3.0;
  const AbaqusCDPSubstepIntegrator deferred(
      state_integrator, {16, strain_limit, 1.0e-8, true});
  const AbaqusCDPSubstepIntegrator legacy(
      state_integrator, {16, strain_limit, 1.0e-8, false});

  const auto result = deferred.integrate({}, target, 1.0e-3, {});
  EXPECT_EQ(result.accepted_substeps, 4u);

  AbaqusCDPLocalIntegrator::State manual_backbone_state;
  std::optional<AbaqusCDPLocalIntegrator::Result> manual_backbone;
  for (unsigned int i = 1; i <= 4; ++i)
  {
    manual_backbone = state_integrator.integrateBackbone(
        substepScale(target, static_cast<double>(i) / 4.0), manual_backbone_state);
    manual_backbone_state = manual_backbone->state;
  }
  ASSERT_TRUE(manual_backbone.has_value());
  const auto expected =
      state_integrator.assembleBackboneResult(target, 1.0e-3, {}, *manual_backbone);

  EXPECT_EQ(result.final_result.state.backbone.plastic_strain,
            expected.state.backbone.plastic_strain);
  EXPECT_EQ(result.final_result.state.viscous_plastic_strain,
            expected.state.viscous_plastic_strain);
  EXPECT_DOUBLE_EQ(result.final_result.state.viscous_tension_damage,
                   expected.state.viscous_tension_damage);
  EXPECT_EQ(result.final_result.cauchy_stress, expected.cauchy_stress);

  const auto legacy_result = legacy.integrate({}, target, 1.0e-3, {});
  EXPECT_NE(result.final_result.state.viscous_plastic_strain,
            legacy_result.final_result.state.viscous_plastic_strain);
}

TEST(AbaqusCDPSubstepIntegrator, DeferredViscousAlgorithmicTangentMatchesDirection)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, substepStateParameters(5.0e-4));
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto target = substepUniaxialElasticStrain(1.05 * initial_tension);
  const AbaqusCDPSubstepIntegrator integrator(
      state_integrator, {16, std::abs(target[0]) / 3.0, 1.0e-8, true});
  const AbaqusCDPSubstepIntegrator::SymmetricTensor direction =
      {1.0, -0.2, -0.2, 0.1, 0.0, 0.0};

  const auto algorithmic = integrator.integrateLinearized({}, target, 1.0e-3, {});
  const auto predicted =
      AbaqusCDPSubstepIntegrator::applyTangent(algorithmic.algorithmic_tangent, direction);
  const auto measured =
      integrator.directionalDerivative({}, target, 1.0e-3, {}, direction, 1.0e-8);
  EXPECT_EQ(algorithmic.result.accepted_substeps, 4u);
  for (std::size_t i = 0; i < 6; ++i)
    EXPECT_NEAR(predicted[i], measured[i], 1.0e-2 * std::max(1.0, std::abs(measured[i])));
}

TEST(AbaqusCDPSubstepIntegrator, AggregatedBackboneHistoryUsesAcceptedStepPlasticIncrement)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, substepStateParameters(5.0e-4));
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto target = substepUniaxialElasticStrain(1.05 * initial_tension);
  const double strain_limit = std::abs(target[0]) / 3.0;
  const AbaqusCDPSubstepIntegrator aggregated(
      state_integrator, {16, strain_limit, 1.0e-8, true, true});

  const auto result = aggregated.integrate({}, target, 1.0e-3, {});
  EXPECT_EQ(result.accepted_substeps, 4u);

  AbaqusCDPLocalIntegrator::State working_backbone;
  std::optional<AbaqusCDPLocalIntegrator::Result> final_backbone;
  for (unsigned int i = 1; i <= 4; ++i)
  {
    final_backbone = state_integrator.integrateBackbone(
        substepScale(target, static_cast<double>(i) / 4.0), working_backbone);
    working_backbone = final_backbone->state;
  }
  ASSERT_TRUE(final_backbone.has_value());
  const auto expected_backbone =
      state_integrator.aggregateBackboneHistory({}, *final_backbone);
  const auto expected =
      state_integrator.assembleBackboneResult(target, 1.0e-3, {}, expected_backbone);

  EXPECT_EQ(result.final_result.state.backbone.plastic_strain,
            expected.state.backbone.plastic_strain);
  EXPECT_DOUBLE_EQ(result.final_result.state.backbone.tensile_equivalent_plastic_strain,
                   expected.state.backbone.tensile_equivalent_plastic_strain);
  EXPECT_DOUBLE_EQ(result.final_result.state.backbone.compressive_equivalent_plastic_strain,
                   expected.state.backbone.compressive_equivalent_plastic_strain);
  EXPECT_DOUBLE_EQ(result.final_result.state.viscous_tension_damage,
                   expected.state.viscous_tension_damage);
  EXPECT_EQ(result.final_result.cauchy_stress, expected.cauchy_stress);
}

TEST(AbaqusCDPSubstepIntegrator, AggregatedBackboneHistoryReferenceTangentMatchesDirection)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, substepStateParameters(5.0e-4));
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto target = substepUniaxialElasticStrain(1.05 * initial_tension);
  const AbaqusCDPSubstepIntegrator integrator(
      state_integrator, {16, std::abs(target[0]) / 3.0, 1.0e-8, true, true});
  const AbaqusCDPSubstepIntegrator::SymmetricTensor direction =
      {1.0, -0.2, -0.2, 0.1, 0.0, 0.0};

  const auto algorithmic = integrator.integrateLinearized({}, target, 1.0e-3, {});
  const auto predicted =
      AbaqusCDPSubstepIntegrator::applyTangent(algorithmic.algorithmic_tangent, direction);
  const auto measured =
      integrator.directionalDerivative({}, target, 1.0e-3, {}, direction, 1.0e-8);
  EXPECT_EQ(algorithmic.result.accepted_substeps, 4u);
  for (std::size_t i = 0; i < 6; ++i)
    EXPECT_NEAR(predicted[i], measured[i], 1.0e-2 * std::max(1.0, std::abs(measured[i])));
}

TEST(AbaqusCDPSubstepIntegrator, AcceptedStepHistoryUsesPlasticWeightedSubstepPathWeights)
{
  const auto table = capturedHistoryTable();
  const AbaqusCDPLocalIntegrator local(
      table, {2.97915e10, 0.2, 36.0, 0.1, 1.16, 0.667, 40, 1.0e-9, 1.0e-7, 1.0e-6});
  const AbaqusCDPStateIntegrator state_integrator(local, {0.0, 1.0, 5.0e-4, 1.0e-12});
  const AbaqusCDPSubstepIntegrator::SymmetricTensor old_strain = {
      5.3299999999999995e-05,
      5.4550000000000005e-05,
      -2.5200000000000000e-04,
      6.8250000000002650e-09,
      -1.7362500000000293e-07,
      2.6949999999999980e-07};
  const AbaqusCDPSubstepIntegrator::SymmetricTensor target_strain = {
      7.360602500000001e-05,
      7.517099000000002e-05,
      -3.366179300000000e-04,
      2.3822382499999515e-08,
      -3.804546000000022e-07,
      5.034524500000326e-08};
  const auto total_increment = substepDifference(target_strain, old_strain);
  const double strain_limit =
      *std::max_element(total_increment.begin(), total_increment.end(), [](double left, double right) {
        return std::abs(left) < std::abs(right);
      }) /
      3.0;
  const AbaqusCDPSubstepIntegrator path_integrated(
      state_integrator, {16, std::abs(strain_limit), 1.0e-8, true, true, true});

  const auto result = path_integrated.integrate(old_strain, target_strain, 1.0e-2, {});
  EXPECT_EQ(result.accepted_substeps, 4u);

  AbaqusCDPLocalIntegrator::State working;
  std::optional<AbaqusCDPLocalIntegrator::Result> final_backbone;
  AbaqusCDPSubstepIntegrator::SymmetricTensor previous_effective_stress =
      state_integrator.backboneEffectiveStress(old_strain, {});
  double tensile_measure = 0.0;
  double tensile_weighted_measure = 0.0;
  double compressive_measure = 0.0;
  double compressive_weighted_measure = 0.0;
  for (unsigned int i = 1; i <= 4; ++i)
  {
    auto target = old_strain;
    for (std::size_t component = 0; component < target.size(); ++component)
      target[component] += static_cast<double>(i) / 4.0 * total_increment[component];
    const auto previous = working;
    final_backbone = state_integrator.integrateBackbone(target, working);
    const auto plastic_increment =
        substepDifference(final_backbone->state.plastic_strain, previous.plastic_strain);
    const auto plastic_principal =
        AbaqusCDPFormula::stressInvariants(plastic_increment).principal_stress;
    const double positive_measure = std::max(0.0, plastic_principal.back());
    const double negative_measure = std::max(0.0, -plastic_principal.front());
    const double start_weight =
        AbaqusCDPFormula::stressInvariants(previous_effective_stress).tension_weight;
    const double end_weight =
        AbaqusCDPFormula::stressInvariants(final_backbone->effective_stress).tension_weight;
    const double average_weight = 0.5 * (start_weight + end_weight);
    tensile_measure += positive_measure;
    tensile_weighted_measure += average_weight * positive_measure;
    compressive_measure += negative_measure;
    compressive_weighted_measure += (1.0 - average_weight) * negative_measure;
    previous_effective_stress = final_backbone->effective_stress;
    working = final_backbone->state;
  }
  ASSERT_TRUE(final_backbone.has_value());
  ASSERT_GT(tensile_measure, 0.0);
  ASSERT_GT(compressive_measure, 0.0);
  const auto expected_backbone = state_integrator.aggregateBackboneHistory(
      {},
      *final_backbone,
      tensile_weighted_measure / tensile_measure,
      compressive_weighted_measure / compressive_measure);
  const auto expected =
      state_integrator.assembleBackboneResult(target_strain, 1.0e-2, {}, expected_backbone);

  EXPECT_DOUBLE_EQ(result.final_result.state.backbone.tensile_equivalent_plastic_strain,
                   expected.state.backbone.tensile_equivalent_plastic_strain);
  EXPECT_DOUBLE_EQ(result.final_result.state.backbone.compressive_equivalent_plastic_strain,
                   expected.state.backbone.compressive_equivalent_plastic_strain);
  EXPECT_EQ(result.final_result.cauchy_stress, expected.cauchy_stress);
}

TEST(AbaqusCDPSubstepIntegrator, AggregatedBackboneHistoryRequiresDeferredViscousState)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, substepStateParameters(5.0e-4));
  EXPECT_THROW(AbaqusCDPSubstepIntegrator(state_integrator, {16, 0.0, 1.0e-8, false, true}),
               std::runtime_error);
  EXPECT_THROW(
      AbaqusCDPSubstepIntegrator(state_integrator, {16, 0.0, 1.0e-8, true, false, true}),
      std::runtime_error);
}

TEST(AbaqusCDPSubstepIntegrator, ExhaustedCutbacksLeaveCallerStateUntouched)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator failing_local(table, substepLocalParameters(0));
  const AbaqusCDPStateIntegrator state_integrator(failing_local, substepStateParameters());
  const AbaqusCDPSubstepIntegrator integrator(state_integrator, {4, 0.0, 1.0e-8});
  AbaqusCDPSubstepIntegrator::State old_state;
  const auto snapshot = old_state;
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;

  try
  {
    integrator.integrate(
        {}, substepUniaxialElasticStrain(1.05 * initial_tension), 1.0e-3, old_state);
    FAIL() << "expected exhausted cutbacks";
  }
  catch (const std::runtime_error & error)
  {
    const std::string diagnostic = error.what();
    EXPECT_NE(diagnostic.find("partition=4"), std::string::npos);
    EXPECT_NE(diagnostic.find("failed_substep="), std::string::npos);
    EXPECT_NE(diagnostic.find("last_error="), std::string::npos);
  }
  EXPECT_EQ(old_state.backbone.plastic_strain, snapshot.backbone.plastic_strain);
  EXPECT_EQ(old_state.viscous_plastic_strain, snapshot.viscous_plastic_strain);
}

TEST(AbaqusCDPSubstepIntegrator, ElasticReferenceTangentMatchesIsotropicTensor)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, substepStateParameters());
  const AbaqusCDPSubstepIntegrator integrator(state_integrator, {8, 0.0, 1.0e-8});
  const auto tangent =
      integrator.referenceTangent({}, substepUniaxialElasticStrain(1.0e6), 1.0e-3, {});
  const double shear =
      substep_youngs_modulus / (2.0 * (1.0 + substep_poissons_ratio));
  const double lambda = substep_youngs_modulus * substep_poissons_ratio /
                        ((1.0 + substep_poissons_ratio) *
                         (1.0 - 2.0 * substep_poissons_ratio));
  const double tolerance = 1.0e-8 * substep_youngs_modulus;

  EXPECT_NEAR(tangent.value[0][0], lambda + 2.0 * shear, tolerance);
  EXPECT_NEAR(tangent.value[0][1], lambda, tolerance);
  EXPECT_NEAR(tangent.value[1][0], lambda, tolerance);
  EXPECT_NEAR(tangent.value[3][3], 2.0 * shear, tolerance);
  EXPECT_NEAR(tangent.value[4][4], 2.0 * shear, tolerance);
  EXPECT_NEAR(tangent.value[5][5], 2.0 * shear, tolerance);
}

TEST(AbaqusCDPSubstepIntegrator, PlasticReferenceTangentMatchesIndependentDirection)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, substepStateParameters(5.0e-4));
  const double perturbation = 1.0e-8;
  const AbaqusCDPSubstepIntegrator integrator(state_integrator, {8, 0.0, perturbation});
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto target = substepUniaxialElasticStrain(1.05 * initial_tension);
  const AbaqusCDPSubstepIntegrator::SymmetricTensor direction =
      {1.0, -0.2, -0.2, 0.1, 0.0, 0.0};

  const auto tangent = integrator.referenceTangent({}, target, 1.0e-3, {});
  const auto predicted = AbaqusCDPSubstepIntegrator::applyTangent(tangent.value, direction);
  const auto measured =
      integrator.directionalDerivative({}, target, 1.0e-3, {}, direction, perturbation);
  for (std::size_t i = 0; i < 6; ++i)
    EXPECT_NEAR(predicted[i], measured[i], 5.0e-3 * std::max(1.0, std::abs(measured[i])));
}

TEST(AbaqusCDPSubstepIntegrator, SingleStepAlgorithmicTangentMatchesReference)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, substepStateParameters(5.0e-4));
  const AbaqusCDPSubstepIntegrator integrator(state_integrator, {8, 0.0, 1.0e-8});
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto target = substepUniaxialElasticStrain(1.05 * initial_tension);
  const AbaqusCDPSubstepIntegrator::SymmetricTensor direction =
      {1.0, -0.2, -0.2, 0.1, 0.0, 0.0};

  const auto algorithmic = integrator.integrateLinearized({}, target, 1.0e-3, {});
  const auto predicted =
      AbaqusCDPSubstepIntegrator::applyTangent(algorithmic.algorithmic_tangent, direction);
  const auto measured =
      integrator.directionalDerivative({}, target, 1.0e-3, {}, direction, 1.0e-8);
  EXPECT_EQ(algorithmic.result.accepted_substeps, 1u);
  for (std::size_t i = 0; i < 6; ++i)
    EXPECT_NEAR(predicted[i], measured[i], 1.0e-2 * std::max(1.0, std::abs(measured[i])));
}

TEST(AbaqusCDPSubstepIntegrator, SubstepAlgorithmicTangentPropagatesStateSensitivity)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, substepStateParameters(5.0e-4));
  const double initial_tension =
      table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0).stress.value;
  const auto target = substepUniaxialElasticStrain(1.05 * initial_tension);
  const AbaqusCDPSubstepIntegrator integrator(
      state_integrator, {16, std::abs(target[0]) / 3.0, 1.0e-8});
  const AbaqusCDPSubstepIntegrator::SymmetricTensor direction =
      {1.0, -0.2, -0.2, 0.1, 0.0, 0.0};

  const auto algorithmic = integrator.integrateLinearized({}, target, 1.0e-3, {});
  const auto predicted =
      AbaqusCDPSubstepIntegrator::applyTangent(algorithmic.algorithmic_tangent, direction);
  const auto measured =
      integrator.directionalDerivative({}, target, 1.0e-3, {}, direction, 1.0e-8);
  EXPECT_EQ(algorithmic.result.accepted_substeps, 4u);
  for (std::size_t i = 0; i < 6; ++i)
    EXPECT_NEAR(predicted[i], measured[i], 1.0e-2 * std::max(1.0, std::abs(measured[i])));
}

TEST(AbaqusCDPSubstepIntegrator, RejectsInvalidConfigurationAndDirection)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state_integrator(local, substepStateParameters());
  EXPECT_THROW(AbaqusCDPSubstepIntegrator(state_integrator, {3, 0.0, 1.0e-8}),
               std::runtime_error);
  EXPECT_THROW(AbaqusCDPSubstepIntegrator(state_integrator, {4, -1.0, 1.0e-8}),
               std::runtime_error);

  const AbaqusCDPSubstepIntegrator integrator(state_integrator, {4, 0.0, 1.0e-8});
  EXPECT_THROW(integrator.directionalDerivative(
                   {}, {}, 1.0e-3, {}, {1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 0.0),
               std::runtime_error);
}

TEST(AbaqusCDPSubstepIntegrator, DiagnosticsRetainFailedPartitionsAndPreserveState)
{
  const auto table=substepReferenceTable();
  const AbaqusCDPLocalIntegrator failing_local(table,substepLocalParameters(0));
  const AbaqusCDPStateIntegrator state(failing_local,substepStateParameters());
  const AbaqusCDPSubstepIntegrator integrator(state,{4,0.0,1e-8});
  AbaqusCDPSubstepIntegrator::State old;
  const auto snapshot=old;
  const double ft=table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION,0).stress.value;
  CDPDiagnostics::Counters costs{};
  CDPDiagnostics::Context context;context.counters=&costs;
  {
    CDPDiagnostics::Binding binding(&context);
    EXPECT_THROW(integrator.integrateLinearized({},substepUniaxialElasticStrain(1.05*ft),1e-3,old),std::runtime_error);
  }
  EXPECT_EQ(CDPDiagnostics::current,nullptr);
  EXPECT_EQ(costs[CDPDiagnostics::PARTITION].calls,3u);
  EXPECT_EQ(costs[CDPDiagnostics::PARTITION].failed,3u);
  EXPECT_GE(costs[CDPDiagnostics::LOCAL].failed,3u);
  EXPECT_GT(costs[CDPDiagnostics::LOCAL].calls,costs[CDPDiagnostics::LOCAL].failed);
  EXPECT_EQ(old.backbone.plastic_strain,snapshot.backbone.plastic_strain);
  EXPECT_EQ(old.viscous_plastic_strain,snapshot.viscous_plastic_strain);
}

TEST(AbaqusCDPSubstepIntegrator, DiagnosticsDoNotChangePlasticMapOrTangent)
{
  const auto table=substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table,substepLocalParameters());
  const AbaqusCDPStateIntegrator state(local,substepStateParameters(5e-4));
  const double ft=table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION,0).stress.value;
  const auto target=substepUniaxialElasticStrain(1.05*ft);
  const AbaqusCDPSubstepIntegrator integrator(state,{16,std::abs(target[0])/3,1e-8});
  const auto reference=integrator.integrateLinearized({},target,1e-3,{});
  CDPDiagnostics::Counters costs{};
  std::ostringstream trace;
  CDPDiagnostics::Context context;context.counters=&costs;context.stream=&trace;context.trace=true;
  CDPDiagnostics::Binding binding(&context);
  const auto observed=integrator.integrateLinearized({},target,1e-3,{});
  EXPECT_EQ(reference.algorithmic_tangent,observed.algorithmic_tangent);
  EXPECT_EQ(reference.result.final_result.cauchy_stress,observed.result.final_result.cauchy_stress);
  EXPECT_EQ(reference.result.final_result.state.backbone.plastic_strain,observed.result.final_result.state.backbone.plastic_strain);
  EXPECT_GT(costs[CDPDiagnostics::AD_JACOBIAN].calls,0u);
  EXPECT_GT(costs[CDPDiagnostics::STATE_CHAIN].calls,0u);
  const auto trace_text=trace.str();
  EXPECT_NE(trace_text.find("\"plastic_increment\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"plastic_increment_principal\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"positive_plastic_principal_increment\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"delta_kappa_t\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"delta_kappa_t_end_rule\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"delta_kappa_c\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"history_integration_rule\":\"end\""),std::string::npos);
  EXPECT_NE(trace_text.find("\"start_target\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"start_backbone_stress\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"start_principal_stress\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"end_principal_stress\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"start_backbone_tension_weight\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"backbone_tension_weight\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"end_backbone_tension_weight\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"maximum_principal_zero_crossing\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"maximum_principal_zero_crossing_fraction\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"viscous_tension_weight\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"active_branch\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"backbone_damage_t\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"backbone_damage_c\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"viscous_damage_t\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"viscous_damage_c\":"),std::string::npos);
  const auto before=costs[CDPDiagnostics::LOCAL].calls;
  {
    CDPDiagnostics::Binding excluded(nullptr);
    integrator.integrateLinearized({},target,1e-3,{});
  }
  EXPECT_EQ(costs[CDPDiagnostics::LOCAL].calls,before);
}

TEST(AbaqusCDPSubstepIntegrator, ReproducesCapturedFirstMixedHistoryState)
{
  const auto table = capturedHistoryTable();
  const AbaqusCDPLocalIntegrator local(
      table, {2.97915e10, 0.2, 36.0, 0.1, 1.16, 0.667, 40, 1.0e-9, 1.0e-7, 1.0e-6});
  const AbaqusCDPStateIntegrator state(local, {0.0, 1.0, 5.0e-4, 1.0e-12});
  const AbaqusCDPSubstepIntegrator integrator(state, {256, 2.5e-5, 1.0e-8});
  const AbaqusCDPSubstepIntegrator::SymmetricTensor old_strain = {
      5.3299999999999995e-05,
      5.4550000000000005e-05,
      -2.5200000000000000e-04,
      6.8250000000002650e-09,
      -1.7362500000000293e-07,
      2.6949999999999980e-07};
  const AbaqusCDPSubstepIntegrator::SymmetricTensor target_strain = {
      5.5330000000000010e-05,
      5.6615000000000006e-05,
      -2.6050000000000004e-04,
      8.5225000000008130e-09,
      -1.9431250000000086e-07,
      2.4759999999999684e-07};

  const auto captured = integrator.integrateLinearized(
      old_strain, target_strain, 5.000000000000004e-4, {});
  const auto repeated = integrator.integrateLinearized(
      old_strain, target_strain, 5.000000000000004e-4, {});
  const auto & result = captured.result.final_result;

  EXPECT_EQ(captured.result.accepted_substeps, 1u);
  EXPECT_EQ(result.backbone.active_branch, AbaqusCDPLocalIntegrator::ActiveBranch::MIXED);
  EXPECT_EQ(captured.algorithmic_tangent, repeated.algorithmic_tangent);
  EXPECT_EQ(result.state.backbone.plastic_strain,
            repeated.result.final_result.state.backbone.plastic_strain);
  EXPECT_NEAR(result.state.backbone.tensile_equivalent_plastic_strain,
              1.7178914856632633e-8,
              1.0e-18);
  EXPECT_NEAR(result.state.backbone.compressive_equivalent_plastic_strain,
              4.6402351272531570e-7,
              1.0e-17);
  EXPECT_NEAR(result.backbone.backbone_tension_damage, 1.491809600909809e-4, 1.0e-14);
  EXPECT_NEAR(result.backbone.backbone_compression_damage, 1.0516310501558594e-3, 1.0e-13);
  EXPECT_NEAR(result.state.viscous_tension_damage, 7.459048004549048e-5, 1.0e-14);
  EXPECT_NEAR(result.state.viscous_compression_damage, 5.258155250779299e-4, 1.0e-14);
  EXPECT_NEAR(result.backbone.residual_norm, 1.4034050595512285e-12, 1.0e-12);

  const AbaqusCDPSubstepIntegrator::SymmetricTensor expected_effective_stress = {
      136541.21892040689,
      168394.9119111658,
      -7692526.374286844,
      211.26311168386522,
      -4816.786551958619,
      6137.723256429379};
  for (std::size_t i = 0; i < expected_effective_stress.size(); ++i)
    EXPECT_NEAR(result.viscous_effective_stress[i], expected_effective_stress[i], 1.0e-6);
}
