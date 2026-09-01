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

CDPMaterialTable
substepReferenceTable()
{
  return CDPMaterialTable(substep_data_dir + "compression_hardening.csv",
                          substep_data_dir + "compression_damage.csv",
                          substep_data_dir + "tension_stiffening.csv",
                          substep_data_dir + "tension_damage.csv",
                          substep_youngs_modulus);
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
  EXPECT_NE(trace_text.find("\"delta_kappa_t\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"delta_kappa_c\":"),std::string::npos);
  EXPECT_NE(trace_text.find("\"backbone_tension_weight\":"),std::string::npos);
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
