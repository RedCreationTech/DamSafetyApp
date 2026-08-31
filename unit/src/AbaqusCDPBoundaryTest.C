#include "AbaqusCDPSubstepIntegrator.h"
#include "gtest/gtest.h"
#include "nlohmann/json.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
using Json = nlohmann::json;
using Tensor = AbaqusCDPLocalIntegrator::SymmetricTensor;
const std::string table_dir = "test/tests/cdp_material_table/uniaxial_tension_20260830/";

Json
failedStates()
{
  std::ifstream input("test/tests/cdp_material_table/compression_boundary_20260831/failed_states.json");
  if (!input)
    throw std::runtime_error("Cannot read frozen CDP boundary fixture");
  return Json::parse(input).at("records");
}

CDPMaterialTable
boundaryTable()
{
  return CDPMaterialTable(table_dir + "compression_hardening.csv",
                          table_dir + "compression_damage.csv",
                          table_dir + "tension_stiffening.csv",
                          table_dir + "tension_damage.csv",
                          2.97915e10);
}

AbaqusCDPLocalIntegrator::Parameters
boundaryParameters(const bool automatic = true)
{
  AbaqusCDPLocalIntegrator::Parameters p{2.97915e10, 0.2, 36.0, 0.1, 1.16, 0.667};
  p.use_automatic_differentiation_jacobian = automatic;
  return p;
}

AbaqusCDPLocalIntegrator::State
localOld(const Json & sample)
{
  return {sample.at("old_plastic").get<Tensor>(),
          sample.at("old_kappa_t").get<double>(),
          sample.at("old_kappa_c").get<double>()};
}

AbaqusCDPStateIntegrator::State
materialOld(const Json & sample)
{
  const auto & s = sample.at("old_state");
  return {{s.at("plastic").get<Tensor>(), s.at("kappa_t"), s.at("kappa_c")},
          s.at("viscous_plastic").get<Tensor>(), s.at("damage_t"), s.at("damage_c")};
}
}

TEST(AbaqusCDPBoundary, FrozenLocalFailuresReproduceWithLegacyADAndFD)
{
  const auto table = boundaryTable();
  unsigned int samples = 0;
  for (const auto & sample : failedStates())
    if (sample.at("event") == "local_failure_jacobian")
    {
      SCOPED_TRACE(sample.at("element").dump() + "/" + sample.at("qp").dump());
      const auto old = localOld(sample);
      for (const bool automatic : {true, false})
      {
        const AbaqusCDPLocalIntegrator local(table, boundaryParameters(automatic));
        EXPECT_THROW(local.integrate(sample.at("target").get<Tensor>(), old), std::runtime_error);
        EXPECT_EQ(old.plastic_strain, localOld(sample).plastic_strain);
        EXPECT_EQ(old.tensile_equivalent_plastic_strain, sample.at("old_kappa_t"));
      }
      ++samples;
    }
  EXPECT_EQ(samples, 16u);
}

TEST(AbaqusCDPBoundary, FrozenCompleteMaterialFailuresReproduceWithLegacyPath)
{
  const auto table = boundaryTable();
  const AbaqusCDPLocalIntegrator local(table, boundaryParameters());
  const AbaqusCDPStateIntegrator state(local, {0.0, 1.0, 5e-4});
  const AbaqusCDPSubstepIntegrator substeps(state, {256, 2.5e-5, 1e-8});
  unsigned int samples = 0;
  for (const auto & sample : failedStates())
    if (sample.at("event") == "material_failure")
    {
      SCOPED_TRACE(sample.at("element").dump() + "/" + sample.at("qp").dump());
      const auto old = materialOld(sample);
      EXPECT_THROW(substeps.integrateLinearized(sample.at("old_strain").get<Tensor>(),
                                                sample.at("new_strain").get<Tensor>(),
                                                sample.at("dt"), old), std::runtime_error);
      ++samples;
    }
  EXPECT_EQ(samples, 13u);
}

TEST(AbaqusCDPBoundary, FeasibleTrialsConvergeFrozenLocalStatesWithoutRelaxingResidual)
{
  const auto table = boundaryTable();
  for (const auto & sample : failedStates())
    if (sample.at("event") == "local_failure_jacobian")
    {
      SCOPED_TRACE(sample.at("element").dump() + "/" + sample.at("qp").dump());
      const auto old = localOld(sample);
      for (const bool automatic : {true, false})
      {
        auto p = boundaryParameters(automatic);
        p.use_bound_feasible_line_search = true;
        const AbaqusCDPLocalIntegrator local(table, p);
        const auto target = sample.at("target").get<Tensor>();
        const auto result = local.integrate(target, old);
        EXPECT_LE(result.residual_norm, 1e-9);
        EXPECT_GE(result.plastic_multiplier, 0.0);
        EXPECT_GE(result.state.tensile_equivalent_plastic_strain, old.tensile_equivalent_plastic_strain);
        EXPECT_GE(result.state.compressive_equivalent_plastic_strain, old.compressive_equivalent_plastic_strain);
        const auto repeated = local.integrate(target, old);
        EXPECT_EQ(result.effective_stress, repeated.effective_stress);
        EXPECT_EQ(result.state.plastic_strain, repeated.state.plastic_strain);
        EXPECT_EQ(old.plastic_strain, localOld(sample).plastic_strain);
      }
    }
}

TEST(AbaqusCDPBoundary, FeasibleTrialsAdvanceFrozenMaterialStatesTransactionally)
{
  const auto table = boundaryTable();
  auto p = boundaryParameters();
  p.use_bound_feasible_line_search = true;
  const AbaqusCDPLocalIntegrator local(table, p);
  const AbaqusCDPStateIntegrator state(local, {0.0, 1.0, 5e-4});
  const AbaqusCDPSubstepIntegrator substeps(state, {256, 2.5e-5, 1e-8});
  for (const auto & sample : failedStates())
    if (sample.at("event") == "material_failure")
    {
      SCOPED_TRACE(sample.at("element").dump() + "/" + sample.at("qp").dump());
      const auto old = materialOld(sample);
      const auto start = sample.at("old_strain").get<Tensor>();
      const auto target = sample.at("new_strain").get<Tensor>();
      const double dt = sample.at("dt");
      const auto result = substeps.integrateLinearized(start, target, dt, old);
      const auto repeated = substeps.integrateLinearized(start, target, dt, old);
      EXPECT_LE(result.result.final_result.backbone.residual_norm, 1e-9);
      EXPECT_EQ(result.result.final_result.cauchy_stress, repeated.result.final_result.cauchy_stress);
      EXPECT_EQ(result.algorithmic_tangent, repeated.algorithmic_tangent);
      EXPECT_EQ(old.backbone.plastic_strain, materialOld(sample).backbone.plastic_strain);
      EXPECT_EQ(old.viscous_plastic_strain, materialOld(sample).viscous_plastic_strain);
    }
}

TEST(AbaqusCDPBoundary, FrozenMaterialTangentsMatchFixedStateDifferenceSweep)
{
  const auto table = boundaryTable();
  auto p = boundaryParameters();
  p.use_bound_feasible_line_search = true;
  const AbaqusCDPLocalIntegrator local(table, p);
  const AbaqusCDPStateIntegrator state(local, {0.0, 1.0, 5e-4});
  const AbaqusCDPSubstepIntegrator substeps(state, {256, 2.5e-5, 1e-8});
  for (const auto & sample : failedStates())
    if (sample.at("event") == "material_failure")
    {
      SCOPED_TRACE(sample.at("element").dump() + "/" + sample.at("qp").dump());
      const auto old = materialOld(sample);
      const auto start = sample.at("old_strain").get<Tensor>();
      const auto target = sample.at("new_strain").get<Tensor>();
      const double dt = sample.at("dt");
      const auto base = substeps.integrateLinearized(start, target, dt, old);
      unsigned int passing_scales = 0;
      for (const double h : {1e-8, 1e-9, 1e-10})
      {
        double worst = 0;
        unsigned int matching_directions = 0;
        for (unsigned int column = 0; column < 6; ++column)
        {
          auto plus = target, minus = target;
          plus[column] += h;
          minus[column] -= h;
          const auto positive = substeps.integrateLinearized(start, plus, dt, old);
          const auto negative = substeps.integrateLinearized(start, minus, dt, old);
          if (positive.result.accepted_substeps != base.result.accepted_substeps ||
              negative.result.accepted_substeps != base.result.accepted_substeps ||
              positive.result.final_result.backbone.active_branch != base.result.final_result.backbone.active_branch ||
              negative.result.final_result.backbone.active_branch != base.result.final_result.backbone.active_branch)
            continue;
          ++matching_directions;
          double error2 = 0, reference2 = 0;
          for (unsigned int row = 0; row < 6; ++row)
          {
            const double fd = (positive.result.final_result.cauchy_stress[row] -
                               negative.result.final_result.cauchy_stress[row]) / (2 * h);
            const double difference = fd - base.algorithmic_tangent[column][row];
            error2 += difference * difference;
            reference2 += fd * fd;
          }
          worst = std::max(worst, std::sqrt(error2 / std::max(reference2, 1.0)));
        }
        std::cout << "BOUND_TANGENT element=" << sample.at("element") << " qp=" << sample.at("qp")
                  << " h=" << h << " directions=" << matching_directions << " error=" << worst << '\n';
        if (matching_directions == 6 && worst <= 1e-4)
          ++passing_scales;
      }
      EXPECT_GE(passing_scales, 1u);
    }
}

TEST(AbaqusCDPBoundary, FeasibleTrialFailureRollbackAndRotatedReplay)
{
  const auto table = boundaryTable();
  auto p = boundaryParameters();
  p.use_bound_feasible_line_search = true;
  const AbaqusCDPLocalIntegrator local(table, p);
  const AbaqusCDPStateIntegrator state(local, {0.0, 1.0, 5e-4});
  const AbaqusCDPSubstepIntegrator substeps(state, {256, 2.5e-5, 1e-8});
  p.maximum_iterations = 0;
  const AbaqusCDPLocalIntegrator failing_local(table, p);
  const AbaqusCDPStateIntegrator failing_state(failing_local, {0.0, 1.0, 5e-4});
  const AbaqusCDPSubstepIntegrator failing(failing_state, {256, 2.5e-5, 1e-8});
  const auto rotate = [](const Tensor & t) -> Tensor {
    return {t[1], t[0], t[2], -t[3], t[5], -t[4]};
  };
  for (const auto & sample : failedStates())
    if (sample.at("event") == "material_failure")
    {
      SCOPED_TRACE(sample.at("element").dump() + "/" + sample.at("qp").dump());
      const auto old = materialOld(sample);
      const auto start = sample.at("old_strain").get<Tensor>();
      const auto target = sample.at("new_strain").get<Tensor>();
      const double dt = sample.at("dt");
      const auto base = substeps.integrateLinearized(start, target, dt, old);
      EXPECT_THROW(failing.integrateLinearized(start, target, dt, old), std::runtime_error);
      const auto retry = substeps.integrateLinearized(start, target, dt, old);
      EXPECT_EQ(base.result.final_result.cauchy_stress, retry.result.final_result.cauchy_stress);
      EXPECT_EQ(base.algorithmic_tangent, retry.algorithmic_tangent);
      auto rotated_old = old;
      rotated_old.backbone.plastic_strain = rotate(old.backbone.plastic_strain);
      rotated_old.viscous_plastic_strain = rotate(old.viscous_plastic_strain);
      const auto rotated = substeps.integrateLinearized(rotate(start), rotate(target), dt, rotated_old);
      const auto expected = rotate(base.result.final_result.cauchy_stress);
      double error2 = 0, reference2 = 0;
      for (unsigned int row = 0; row < 6; ++row)
      {
        const double difference = rotated.result.final_result.cauchy_stress[row] - expected[row];
        error2 += difference * difference;
        reference2 += expected[row] * expected[row];
      }
      EXPECT_LE(std::sqrt(error2 / reference2), 1e-6);
      EXPECT_NEAR(rotated.result.final_result.state.viscous_tension_damage,
                  base.result.final_result.state.viscous_tension_damage, 1e-9);
      EXPECT_NEAR(rotated.result.final_result.state.viscous_compression_damage,
                  base.result.final_result.state.viscous_compression_damage, 1e-9);
    }
}
