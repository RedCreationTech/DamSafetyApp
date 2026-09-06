#include "AbaqusCDPLocalIntegrator.h"

#include "gtest/gtest.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <optional>

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

TEST(AbaqusCDPLocalIntegrator, FailureCaptureDoesNotChangeReturnOrFailure)
{
  const auto table = referenceTable();
  auto off = parameters();
  off.minimum_line_search = 1.0; // exercise rejection without a long backtracking path
  auto on = off;
  on.capture_failure_context = true;
  const AbaqusCDPLocalIntegrator plain(table, off), captured(table, on);
  unsigned int payloads = 0;
  for (double strain : {1e-6, 1e-4, 1e-3, 1e-2, -1e-4, -1e-3, -1e-2})
  {
    const AbaqusCDPLocalIntegrator::SymmetricTensor target =
        {strain, -0.2 * strain, -0.2 * strain, 0.0, 0.0, 0.0};
    AbaqusCDPLocalIntegrator::State old_state;
    std::string expected_error;
    std::optional<AbaqusCDPLocalIntegrator::LinearizedResult> reference;
    try
    {
      reference = plain.integrateLinearized(target, old_state);
    }
    catch (const std::runtime_error & error)
    {
      expected_error = error.what();
    }
    if (reference)
    {
      std::optional<AbaqusCDPLocalIntegrator::LinearizedResult> observed;
      ASSERT_NO_THROW(observed = captured.integrateLinearized(target, old_state));
      const auto & a = *reference;
      const auto & b = *observed;
      EXPECT_EQ(a.result.effective_stress, b.result.effective_stress);
      EXPECT_EQ(a.result.state.plastic_strain, b.result.state.plastic_strain);
      EXPECT_EQ(a.result.state.tensile_equivalent_plastic_strain,
                b.result.state.tensile_equivalent_plastic_strain);
      EXPECT_EQ(a.result.state.compressive_equivalent_plastic_strain,
                b.result.state.compressive_equivalent_plastic_strain);
      EXPECT_EQ(a.derivative, b.derivative);
      EXPECT_EQ(a.result.iterations, b.result.iterations);
    }
    if (!expected_error.empty())
    {
      try
      {
        captured.integrateLinearized(target, old_state);
        FAIL() << "capture changed failure into success";
      }
      catch (const std::runtime_error & error)
      {
        const std::string actual = error.what();
        EXPECT_EQ(actual.substr(0, expected_error.size()), expected_error);
        if (expected_error.find("line search failed") != std::string::npos)
        {
          EXPECT_NE(actual.find("CDP_LOCAL_FAILURE_V1"), std::string::npos);
          EXPECT_NE(actual.find("negative_rejections"), std::string::npos);
          EXPECT_NE(actual.find("old_plastic_strain"), std::string::npos);
          ++payloads;
        }
      }
    }
    EXPECT_EQ(old_state.tensile_equivalent_plastic_strain, 0.0);
    EXPECT_EQ(old_state.compressive_equivalent_plastic_strain, 0.0);
  }
  EXPECT_GT(payloads, 0u);
}

TEST(AbaqusCDPLocalIntegrator, CapturedD01Failure0ConvergesWithoutRelaxation)
{
  const std::string d = "test/tests/cdp_material_table/reloading_capture_data/";
  const CDPMaterialTable table(d + "compression_hardening.csv", d + "compression_damage.csv",
                              d + "tension_stiffening.csv", d + "tension_damage.csv", 29791500000.0);
  const AbaqusCDPLocalIntegrator integrator(table, {29791500000.0, 0.2, 36, 0.1, 1.16, 0.667});
  const AbaqusCDPLocalIntegrator::SymmetricTensor target = {-9.0365343635225316e-05, -9.036534363203684e-05, 0.00093614078695955288, 1.2324873920195368e-05, -1.7592418891230148e-05, -1.7592418905879295e-05};
  AbaqusCDPLocalIntegrator::State old;
  old.plastic_strain = {-0.0001700362842915049, -0.00017003628428836188, 0.00084379512253861754, 1.216360837019724e-05, -1.7356925963159499e-05, -1.7356925977610017e-05};
  old.tensile_equivalent_plastic_strain = 0.00084441609852190045;
  old.compressive_equivalent_plastic_strain = 0;
  const auto saved = old;
  const auto linearized = integrator.integrateLinearized(target, old);
  const auto & result = linearized.result;
  EXPECT_LT(result.residual_norm, 1e-9);
  EXPECT_LE(result.iterations, 40u);
  EXPECT_GT(result.plastic_multiplier, 0.0);
  EXPECT_GE(result.state.tensile_equivalent_plastic_strain, old.tensile_equivalent_plastic_strain);
  EXPECT_GE(result.state.compressive_equivalent_plastic_strain, old.compressive_equivalent_plastic_strain);
  EXPECT_EQ(old.plastic_strain, saved.plastic_strain);
  EXPECT_EQ(old.tensile_equivalent_plastic_strain, saved.tensile_equivalent_plastic_strain);
  EXPECT_EQ(old.compressive_equivalent_plastic_strain, saved.compressive_equivalent_plastic_strain);
  auto plus = target, minus = target;
  const double step = 1e-10;
  plus[2] += step;
  minus[2] -= step;
  const auto upper = integrator.integrate(plus, old);
  const auto lower = integrator.integrate(minus, old);
  for (unsigned int k = 0; k < 6; ++k)
  {
    const double finite_difference = (upper.effective_stress[k] - lower.effective_stress[k]) / (2 * step);
    EXPECT_NEAR(linearized.derivative[2][k], finite_difference,
                1e-4 * std::max(1.0e6, std::abs(finite_difference)));
  }
}

TEST(AbaqusCDPLocalIntegrator, CapturedD02Failure0ConvergesWithoutRelaxation)
{
  const std::string d = "test/tests/cdp_material_table/reloading_capture_data/";
  const CDPMaterialTable table(d + "compression_hardening.csv", d + "compression_damage.csv",
                              d + "tension_stiffening.csv", d + "tension_damage.csv", 29791500000.0);
  const AbaqusCDPLocalIntegrator integrator(table, {29791500000.0, 0.2, 36, 0.1, 1.16, 0.667});
  const AbaqusCDPLocalIntegrator::SymmetricTensor target = {0.0015240141114345868, 0.0015240141110805254, -0.0013975892305307224, -0.00021690460335515682, 0.00017610941571070507, 0.0001761094157224077};
  AbaqusCDPLocalIntegrator::State old;
  old.plastic_strain = {0.0014872687374944899, 0.0014872687371291821, -0.0014620339967781621, -0.00020665604500810279, 0.00016646232338713882, 0.00016646232341545116};
  old.tensile_equivalent_plastic_strain = 0.00012298361418062488;
  old.compressive_equivalent_plastic_strain = 0.0014613947559259319;
  const auto saved = old;
  const auto linearized = integrator.integrateLinearized(target, old);
  const auto & result = linearized.result;
  EXPECT_LT(result.residual_norm, 1e-9);
  EXPECT_LE(result.iterations, 40u);
  EXPECT_GT(result.plastic_multiplier, 0.0);
  EXPECT_GE(result.state.tensile_equivalent_plastic_strain, old.tensile_equivalent_plastic_strain);
  EXPECT_GE(result.state.compressive_equivalent_plastic_strain, old.compressive_equivalent_plastic_strain);
  EXPECT_EQ(old.plastic_strain, saved.plastic_strain);
  EXPECT_EQ(old.tensile_equivalent_plastic_strain, saved.tensile_equivalent_plastic_strain);
  EXPECT_EQ(old.compressive_equivalent_plastic_strain, saved.compressive_equivalent_plastic_strain);
  auto plus = target, minus = target;
  const double step = 1e-10;
  plus[2] += step;
  minus[2] -= step;
  const auto upper = integrator.integrate(plus, old);
  const auto lower = integrator.integrate(minus, old);
  for (unsigned int k = 0; k < 6; ++k)
  {
    const double finite_difference = (upper.effective_stress[k] - lower.effective_stress[k]) / (2 * step);
    EXPECT_NEAR(linearized.derivative[2][k], finite_difference,
                1e-4 * std::max(1.0e6, std::abs(finite_difference)));
  }
}

TEST(AbaqusCDPLocalIntegrator, CapturedD02Failure1ConvergesWithoutRelaxation)
{
  const std::string d = "test/tests/cdp_material_table/reloading_capture_data/";
  const CDPMaterialTable table(d + "compression_hardening.csv", d + "compression_damage.csv",
                              d + "tension_stiffening.csv", d + "tension_damage.csv", 29791500000.0);
  const AbaqusCDPLocalIntegrator integrator(table, {29791500000.0, 0.2, 36, 0.1, 1.16, 0.667});
  const AbaqusCDPLocalIntegrator::SymmetricTensor target = {0.0009531850571030951, 0.0009531850548886397, 0.00015815542001929389, 0.00061822768025922, -0.0013701250945158156, -0.0013701250960286499};
  AbaqusCDPLocalIntegrator::State old;
  old.plastic_strain = {0.00089439568582129524, 0.00089439568360668502, 8.5302688470749991e-05, 0.00061818102327210037, -0.0013700165132083122, -0.0013700165147215513};
  old.tensile_equivalent_plastic_strain = 0.0021957686373652988;
  old.compressive_equivalent_plastic_strain = 0.00043273549203351008;
  const auto saved = old;
  const auto linearized = integrator.integrateLinearized(target, old);
  const auto & result = linearized.result;
  EXPECT_LT(result.residual_norm, 1e-9);
  EXPECT_LE(result.iterations, 40u);
  EXPECT_GT(result.plastic_multiplier, 0.0);
  EXPECT_GE(result.state.tensile_equivalent_plastic_strain, old.tensile_equivalent_plastic_strain);
  EXPECT_GE(result.state.compressive_equivalent_plastic_strain, old.compressive_equivalent_plastic_strain);
  EXPECT_EQ(old.plastic_strain, saved.plastic_strain);
  EXPECT_EQ(old.tensile_equivalent_plastic_strain, saved.tensile_equivalent_plastic_strain);
  EXPECT_EQ(old.compressive_equivalent_plastic_strain, saved.compressive_equivalent_plastic_strain);
  auto plus = target, minus = target;
  const double step = 1e-10;
  plus[2] += step;
  minus[2] -= step;
  const auto upper = integrator.integrate(plus, old);
  const auto lower = integrator.integrate(minus, old);
  for (unsigned int k = 0; k < 6; ++k)
  {
    const double finite_difference = (upper.effective_stress[k] - lower.effective_stress[k]) / (2 * step);
    EXPECT_NEAR(linearized.derivative[2][k], finite_difference,
                1e-4 * std::max(1.0e6, std::abs(finite_difference)));
  }
}

TEST(AbaqusCDPLocalIntegrator, CapturedD02Failure2ConvergesWithoutRelaxation)
{
  const std::string d = "test/tests/cdp_material_table/reloading_capture_data/";
  const CDPMaterialTable table(d + "compression_hardening.csv", d + "compression_damage.csv",
                              d + "tension_stiffening.csv", d + "tension_damage.csv", 29791500000.0);
  const AbaqusCDPLocalIntegrator integrator(table, {29791500000.0, 0.2, 36, 0.1, 1.16, 0.667});
  const AbaqusCDPLocalIntegrator::SymmetricTensor target = {0.00096133133028604694, 0.00096133131583884689, 0.00018357006345837883, 0.00062554901815871085, -0.0013831792766431663, -0.0013831792782876442};
  AbaqusCDPLocalIntegrator::State old;
  old.plastic_strain = {0.00089792970657001636, 0.00089792969212419726, 0.00011291180885880218, 0.00062554735417117379, -0.0013831757902785732, -0.0013831757919230423};
  old.tensile_equivalent_plastic_strain = 0.0022492747425780671;
  old.compressive_equivalent_plastic_strain = 0.00043273549203351008;
  const auto saved = old;
  const auto linearized = integrator.integrateLinearized(target, old);
  const auto & result = linearized.result;
  EXPECT_LT(result.residual_norm, 1e-9);
  EXPECT_LE(result.iterations, 40u);
  EXPECT_GT(result.plastic_multiplier, 0.0);
  EXPECT_GE(result.state.tensile_equivalent_plastic_strain, old.tensile_equivalent_plastic_strain);
  EXPECT_GE(result.state.compressive_equivalent_plastic_strain, old.compressive_equivalent_plastic_strain);
  EXPECT_EQ(old.plastic_strain, saved.plastic_strain);
  EXPECT_EQ(old.tensile_equivalent_plastic_strain, saved.tensile_equivalent_plastic_strain);
  EXPECT_EQ(old.compressive_equivalent_plastic_strain, saved.compressive_equivalent_plastic_strain);
  auto plus = target, minus = target;
  const double step = 1e-10;
  plus[2] += step;
  minus[2] -= step;
  const auto upper = integrator.integrate(plus, old);
  const auto lower = integrator.integrate(minus, old);
  for (unsigned int k = 0; k < 6; ++k)
  {
    const double finite_difference = (upper.effective_stress[k] - lower.effective_stress[k]) / (2 * step);
    EXPECT_NEAR(linearized.derivative[2][k], finite_difference,
                1e-4 * std::max(1.0e6, std::abs(finite_difference)));
  }
}

TEST(AbaqusCDPLocalIntegrator, CapturedD02Failure3ConvergesWithoutRelaxation)
{
  const std::string d = "test/tests/cdp_material_table/reloading_capture_data/";
  const CDPMaterialTable table(d + "compression_hardening.csv", d + "compression_damage.csv",
                              d + "tension_stiffening.csv", d + "tension_damage.csv", 29791500000.0);
  const AbaqusCDPLocalIntegrator integrator(table, {29791500000.0, 0.2, 36, 0.1, 1.16, 0.667});
  const AbaqusCDPLocalIntegrator::SymmetricTensor target = {0.000972388249845218, 0.00097238851354146155, 0.00019442585526969455, 0.0006369081687145761, -0.0014009130930480956, -0.0014009128555236452};
  AbaqusCDPLocalIntegrator::State old;
  old.plastic_strain = {0.00090746960819255583, 0.00090746987188716798, 0.00012445839532105251, 0.00063690805969194729, -0.0014009128986647683, -0.0014009126611421736};
  old.tensile_equivalent_plastic_strain = 0.0022945280968675206;
  old.compressive_equivalent_plastic_strain = 0.00043273549203351008;
  const auto saved = old;
  const auto linearized = integrator.integrateLinearized(target, old);
  const auto & result = linearized.result;
  EXPECT_LT(result.residual_norm, 1e-9);
  EXPECT_LE(result.iterations, 40u);
  EXPECT_GT(result.plastic_multiplier, 0.0);
  EXPECT_GE(result.state.tensile_equivalent_plastic_strain, old.tensile_equivalent_plastic_strain);
  EXPECT_GE(result.state.compressive_equivalent_plastic_strain, old.compressive_equivalent_plastic_strain);
  EXPECT_EQ(old.plastic_strain, saved.plastic_strain);
  EXPECT_EQ(old.tensile_equivalent_plastic_strain, saved.tensile_equivalent_plastic_strain);
  EXPECT_EQ(old.compressive_equivalent_plastic_strain, saved.compressive_equivalent_plastic_strain);
  auto plus = target, minus = target;
  const double step = 1e-10;
  plus[2] += step;
  minus[2] -= step;
  const auto upper = integrator.integrate(plus, old);
  const auto lower = integrator.integrate(minus, old);
  for (unsigned int k = 0; k < 6; ++k)
  {
    const double finite_difference = (upper.effective_stress[k] - lower.effective_stress[k]) / (2 * step);
    EXPECT_NEAR(linearized.derivative[2][k], finite_difference,
                1e-4 * std::max(1.0e6, std::abs(finite_difference)));
  }
}
