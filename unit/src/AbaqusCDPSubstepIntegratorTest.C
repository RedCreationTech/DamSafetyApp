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

TEST(AbaqusCDPSubstepIntegrator, TerminalRetryLeavesOriginalSuccessfulPartitionsAndTangentsExact)
{
  const auto table = substepReferenceTable();
  const AbaqusCDPLocalIntegrator local(table, substepLocalParameters());
  const AbaqusCDPStateIntegrator state(local, substepStateParameters(5e-4));
  for (double stress : {2.1e6, -2.5e7})
    for (unsigned int count : {1u, 2u, 4u, 8u, 16u})
    {
      const auto target = substepUniaxialElasticStrain(stress);
      const double increment = std::abs(target[0]) / count * 1.001;
      const AbaqusCDPSubstepIntegrator before(state, {16, increment, 1e-8, false});
      const AbaqusCDPSubstepIntegrator after(state, {16, increment, 1e-8, true});
      const auto a = before.integrateLinearized({}, target, 1e-3, {});
      const auto b = after.integrateLinearized({}, target, 1e-3, {});
      EXPECT_EQ(a.result.accepted_substeps, b.result.accepted_substeps);
      EXPECT_EQ(a.result.total_local_iterations, b.result.total_local_iterations);
      EXPECT_EQ(a.algorithmic_tangent, b.algorithmic_tangent);
      EXPECT_EQ(a.result.final_result.cauchy_stress, b.result.final_result.cauchy_stress);
      EXPECT_EQ(a.result.final_result.state.backbone.plastic_strain,
                b.result.final_result.state.backbone.plastic_strain);
      EXPECT_EQ(a.result.final_result.state.viscous_plastic_strain,
                b.result.final_result.state.viscous_plastic_strain);
      EXPECT_EQ(a.result.final_result.state.viscous_tension_damage,
                b.result.final_result.state.viscous_tension_damage);
      EXPECT_EQ(a.result.final_result.state.viscous_compression_damage,
                b.result.final_result.state.viscous_compression_damage);
      AbaqusCDPStateIntegrator::State manual;
      std::optional<AbaqusCDPStateIntegrator::Result> last;
      for (unsigned int i = 1; i <= a.result.accepted_substeps; ++i)
      {
        last = state.integrate(substepScale(target, double(i) / a.result.accepted_substeps),
                               1e-3 / a.result.accepted_substeps, manual);
        manual = last->state;
      }
      EXPECT_EQ(last->cauchy_stress, b.result.final_result.cauchy_stress);
      EXPECT_EQ(manual.viscous_plastic_strain, b.result.final_result.state.viscous_plastic_strain);
      EXPECT_EQ(manual.viscous_tension_damage, b.result.final_result.state.viscous_tension_damage);
      EXPECT_EQ(manual.viscous_compression_damage, b.result.final_result.state.viscous_compression_damage);
    }
}

TEST(AbaqusCDPSubstepIntegrator, TerminalRescueD01Capture0IsTransactional)
{
  const std::string d = "test/tests/cdp_material_table/reloading_capture_data/";
  const CDPMaterialTable table(d+"compression_hardening.csv",d+"compression_damage.csv",
                              d+"tension_stiffening.csv",d+"tension_damage.csv",29791500000.0);
  const AbaqusCDPLocalIntegrator local(table, {29791500000.0,0.2,36,0.1,1.16,0.667});
  const AbaqusCDPStateIntegrator state(local, {1,0,5e-4,1e-12});
  AbaqusCDPStateIntegrator::State old;
  old.backbone.plastic_strain={-0.0001700362842915049,-0.00017003628428836188,0.00084379512253861754,1.216360837019724e-05,-1.7356925963159499e-05,-1.7356925977610017e-05};
  old.backbone.tensile_equivalent_plastic_strain=0.00084441609852190045;
  old.backbone.compressive_equivalent_plastic_strain=0;
  old.viscous_plastic_strain=old.backbone.plastic_strain;
  const auto saved=old;
  const AbaqusCDPSubstepIntegrator::SymmetricTensor target={-9.0365343635225316e-05,-9.036534363203684e-05,0.00093614078695955288,1.2324873920195368e-05,-1.7592418891230148e-05,-1.7592418905879295e-05};
  // One permitted partition isolates dispatch of the terminal fallback. The
  // viscous state is a declared software fixture, not a captured accepted IP.
  const AbaqusCDPSubstepIntegrator before(state, {1,0,1e-8,false});
  const AbaqusCDPSubstepIntegrator after(state, {1,0,1e-8,true});
  EXPECT_THROW(before.integrateLinearized(target,target,1e-5,old),std::runtime_error);
  EXPECT_EQ(old.backbone.plastic_strain,saved.backbone.plastic_strain);
  EXPECT_EQ(old.viscous_plastic_strain,saved.viscous_plastic_strain);
  const auto result=after.integrateLinearized(target,target,1e-5,old);
  const auto plain=after.integrate(target,target,1e-5,old);
  EXPECT_EQ(result.result.accepted_substeps,1u);
  EXPECT_LT(result.result.final_result.backbone.residual_norm,1e-9);
  EXPECT_EQ(result.result.final_result.cauchy_stress,plain.final_result.cauchy_stress);
  EXPECT_EQ(old.backbone.tensile_equivalent_plastic_strain,saved.backbone.tensile_equivalent_plastic_strain);
  EXPECT_EQ(old.backbone.compressive_equivalent_plastic_strain,saved.backbone.compressive_equivalent_plastic_strain);
  EXPECT_EQ(old.viscous_tension_damage,saved.viscous_tension_damage);
  EXPECT_EQ(old.viscous_compression_damage,saved.viscous_compression_damage);
}

TEST(AbaqusCDPSubstepIntegrator, TerminalRescueD02Capture0IsTransactional)
{
  const std::string d = "test/tests/cdp_material_table/reloading_capture_data/";
  const CDPMaterialTable table(d+"compression_hardening.csv",d+"compression_damage.csv",
                              d+"tension_stiffening.csv",d+"tension_damage.csv",29791500000.0);
  const AbaqusCDPLocalIntegrator local(table, {29791500000.0,0.2,36,0.1,1.16,0.667});
  const AbaqusCDPStateIntegrator state(local, {1,0,5e-4,1e-12});
  AbaqusCDPStateIntegrator::State old;
  old.backbone.plastic_strain={0.0014872687374944899,0.0014872687371291821,-0.0014620339967781621,-0.00020665604500810279,0.00016646232338713882,0.00016646232341545116};
  old.backbone.tensile_equivalent_plastic_strain=0.00012298361418062488;
  old.backbone.compressive_equivalent_plastic_strain=0.0014613947559259319;
  old.viscous_plastic_strain=old.backbone.plastic_strain;
  const auto saved=old;
  const AbaqusCDPSubstepIntegrator::SymmetricTensor target={0.0015240141114345868,0.0015240141110805254,-0.0013975892305307224,-0.00021690460335515682,0.00017610941571070507,0.0001761094157224077};
  // One permitted partition isolates dispatch of the terminal fallback. The
  // viscous state is a declared software fixture, not a captured accepted IP.
  const AbaqusCDPSubstepIntegrator before(state, {1,0,1e-8,false});
  const AbaqusCDPSubstepIntegrator after(state, {1,0,1e-8,true});
  EXPECT_THROW(before.integrateLinearized(target,target,1e-5,old),std::runtime_error);
  EXPECT_EQ(old.backbone.plastic_strain,saved.backbone.plastic_strain);
  EXPECT_EQ(old.viscous_plastic_strain,saved.viscous_plastic_strain);
  const auto result=after.integrateLinearized(target,target,1e-5,old);
  const auto plain=after.integrate(target,target,1e-5,old);
  EXPECT_EQ(result.result.accepted_substeps,1u);
  EXPECT_LT(result.result.final_result.backbone.residual_norm,1e-9);
  EXPECT_EQ(result.result.final_result.cauchy_stress,plain.final_result.cauchy_stress);
  EXPECT_EQ(old.backbone.tensile_equivalent_plastic_strain,saved.backbone.tensile_equivalent_plastic_strain);
  EXPECT_EQ(old.backbone.compressive_equivalent_plastic_strain,saved.backbone.compressive_equivalent_plastic_strain);
  EXPECT_EQ(old.viscous_tension_damage,saved.viscous_tension_damage);
  EXPECT_EQ(old.viscous_compression_damage,saved.viscous_compression_damage);
}

TEST(AbaqusCDPSubstepIntegrator, TerminalRescueD02Capture1IsTransactional)
{
  const std::string d = "test/tests/cdp_material_table/reloading_capture_data/";
  const CDPMaterialTable table(d+"compression_hardening.csv",d+"compression_damage.csv",
                              d+"tension_stiffening.csv",d+"tension_damage.csv",29791500000.0);
  const AbaqusCDPLocalIntegrator local(table, {29791500000.0,0.2,36,0.1,1.16,0.667});
  const AbaqusCDPStateIntegrator state(local, {1,0,5e-4,1e-12});
  AbaqusCDPStateIntegrator::State old;
  old.backbone.plastic_strain={0.00089439568582129524,0.00089439568360668502,8.5302688470749991e-05,0.00061818102327210037,-0.0013700165132083122,-0.0013700165147215513};
  old.backbone.tensile_equivalent_plastic_strain=0.0021957686373652988;
  old.backbone.compressive_equivalent_plastic_strain=0.00043273549203351008;
  old.viscous_plastic_strain=old.backbone.plastic_strain;
  const auto saved=old;
  const AbaqusCDPSubstepIntegrator::SymmetricTensor target={0.0009531850571030951,0.0009531850548886397,0.00015815542001929389,0.00061822768025922,-0.0013701250945158156,-0.0013701250960286499};
  // One permitted partition isolates dispatch of the terminal fallback. The
  // viscous state is a declared software fixture, not a captured accepted IP.
  const AbaqusCDPSubstepIntegrator before(state, {1,0,1e-8,false});
  const AbaqusCDPSubstepIntegrator after(state, {1,0,1e-8,true});
  EXPECT_THROW(before.integrateLinearized(target,target,1e-5,old),std::runtime_error);
  EXPECT_EQ(old.backbone.plastic_strain,saved.backbone.plastic_strain);
  EXPECT_EQ(old.viscous_plastic_strain,saved.viscous_plastic_strain);
  const auto result=after.integrateLinearized(target,target,1e-5,old);
  const auto plain=after.integrate(target,target,1e-5,old);
  EXPECT_EQ(result.result.accepted_substeps,1u);
  EXPECT_LT(result.result.final_result.backbone.residual_norm,1e-9);
  EXPECT_EQ(result.result.final_result.cauchy_stress,plain.final_result.cauchy_stress);
  EXPECT_EQ(old.backbone.tensile_equivalent_plastic_strain,saved.backbone.tensile_equivalent_plastic_strain);
  EXPECT_EQ(old.backbone.compressive_equivalent_plastic_strain,saved.backbone.compressive_equivalent_plastic_strain);
  EXPECT_EQ(old.viscous_tension_damage,saved.viscous_tension_damage);
  EXPECT_EQ(old.viscous_compression_damage,saved.viscous_compression_damage);
}

TEST(AbaqusCDPSubstepIntegrator, TerminalRescueD02Capture2IsTransactional)
{
  const std::string d = "test/tests/cdp_material_table/reloading_capture_data/";
  const CDPMaterialTable table(d+"compression_hardening.csv",d+"compression_damage.csv",
                              d+"tension_stiffening.csv",d+"tension_damage.csv",29791500000.0);
  const AbaqusCDPLocalIntegrator local(table, {29791500000.0,0.2,36,0.1,1.16,0.667});
  const AbaqusCDPStateIntegrator state(local, {1,0,5e-4,1e-12});
  AbaqusCDPStateIntegrator::State old;
  old.backbone.plastic_strain={0.00089792970657001636,0.00089792969212419726,0.00011291180885880218,0.00062554735417117379,-0.0013831757902785732,-0.0013831757919230423};
  old.backbone.tensile_equivalent_plastic_strain=0.0022492747425780671;
  old.backbone.compressive_equivalent_plastic_strain=0.00043273549203351008;
  old.viscous_plastic_strain=old.backbone.plastic_strain;
  const auto saved=old;
  const AbaqusCDPSubstepIntegrator::SymmetricTensor target={0.00096133133028604694,0.00096133131583884689,0.00018357006345837883,0.00062554901815871085,-0.0013831792766431663,-0.0013831792782876442};
  // One permitted partition isolates dispatch of the terminal fallback. The
  // viscous state is a declared software fixture, not a captured accepted IP.
  const AbaqusCDPSubstepIntegrator before(state, {1,0,1e-8,false});
  const AbaqusCDPSubstepIntegrator after(state, {1,0,1e-8,true});
  EXPECT_THROW(before.integrateLinearized(target,target,1e-5,old),std::runtime_error);
  EXPECT_EQ(old.backbone.plastic_strain,saved.backbone.plastic_strain);
  EXPECT_EQ(old.viscous_plastic_strain,saved.viscous_plastic_strain);
  const auto result=after.integrateLinearized(target,target,1e-5,old);
  const auto plain=after.integrate(target,target,1e-5,old);
  EXPECT_EQ(result.result.accepted_substeps,1u);
  EXPECT_LT(result.result.final_result.backbone.residual_norm,1e-9);
  EXPECT_EQ(result.result.final_result.cauchy_stress,plain.final_result.cauchy_stress);
  EXPECT_EQ(old.backbone.tensile_equivalent_plastic_strain,saved.backbone.tensile_equivalent_plastic_strain);
  EXPECT_EQ(old.backbone.compressive_equivalent_plastic_strain,saved.backbone.compressive_equivalent_plastic_strain);
  EXPECT_EQ(old.viscous_tension_damage,saved.viscous_tension_damage);
  EXPECT_EQ(old.viscous_compression_damage,saved.viscous_compression_damage);
}

TEST(AbaqusCDPSubstepIntegrator, TerminalRescueD02Capture3IsTransactional)
{
  const std::string d = "test/tests/cdp_material_table/reloading_capture_data/";
  const CDPMaterialTable table(d+"compression_hardening.csv",d+"compression_damage.csv",
                              d+"tension_stiffening.csv",d+"tension_damage.csv",29791500000.0);
  const AbaqusCDPLocalIntegrator local(table, {29791500000.0,0.2,36,0.1,1.16,0.667});
  const AbaqusCDPStateIntegrator state(local, {1,0,5e-4,1e-12});
  AbaqusCDPStateIntegrator::State old;
  old.backbone.plastic_strain={0.00090746960819255583,0.00090746987188716798,0.00012445839532105251,0.00063690805969194729,-0.0014009128986647683,-0.0014009126611421736};
  old.backbone.tensile_equivalent_plastic_strain=0.0022945280968675206;
  old.backbone.compressive_equivalent_plastic_strain=0.00043273549203351008;
  old.viscous_plastic_strain=old.backbone.plastic_strain;
  const auto saved=old;
  const AbaqusCDPSubstepIntegrator::SymmetricTensor target={0.000972388249845218,0.00097238851354146155,0.00019442585526969455,0.0006369081687145761,-0.0014009130930480956,-0.0014009128555236452};
  // One permitted partition isolates dispatch of the terminal fallback. The
  // viscous state is a declared software fixture, not a captured accepted IP.
  const AbaqusCDPSubstepIntegrator before(state, {1,0,1e-8,false});
  const AbaqusCDPSubstepIntegrator after(state, {1,0,1e-8,true});
  EXPECT_THROW(before.integrateLinearized(target,target,1e-5,old),std::runtime_error);
  EXPECT_EQ(old.backbone.plastic_strain,saved.backbone.plastic_strain);
  EXPECT_EQ(old.viscous_plastic_strain,saved.viscous_plastic_strain);
  const auto result=after.integrateLinearized(target,target,1e-5,old);
  const auto plain=after.integrate(target,target,1e-5,old);
  EXPECT_EQ(result.result.accepted_substeps,1u);
  EXPECT_LT(result.result.final_result.backbone.residual_norm,1e-9);
  EXPECT_EQ(result.result.final_result.cauchy_stress,plain.final_result.cauchy_stress);
  EXPECT_EQ(old.backbone.tensile_equivalent_plastic_strain,saved.backbone.tensile_equivalent_plastic_strain);
  EXPECT_EQ(old.backbone.compressive_equivalent_plastic_strain,saved.backbone.compressive_equivalent_plastic_strain);
  EXPECT_EQ(old.viscous_tension_damage,saved.viscous_tension_damage);
  EXPECT_EQ(old.viscous_compression_damage,saved.viscous_compression_damage);
}
