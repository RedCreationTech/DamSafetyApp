#include "AbaqusCDPSubstepIntegrator.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace
{
using Six = AbaqusCDPSubstepIntegrator::SymmetricTensor;
using State = AbaqusCDPStateIntegrator::State;
using StateValues = std::array<double, 16>;

constexpr double youngs_modulus = 2.97915e10;
const std::string data_dir = "unit/data/ds04g-sp1/";
const std::array<const char *, 6> tensor_names = {"xx", "yy", "zz", "xy", "yz", "xz"};
const std::array<const char *, 16> state_names = {"backbone_xx",
                                                  "backbone_yy",
                                                  "backbone_zz",
                                                  "backbone_xy",
                                                  "backbone_yz",
                                                  "backbone_xz",
                                                  "kappa_t",
                                                  "kappa_c",
                                                  "viscous_xx",
                                                  "viscous_yy",
                                                  "viscous_zz",
                                                  "viscous_xy",
                                                  "viscous_yz",
                                                  "viscous_xz",
                                                  "damage_t",
                                                  "damage_c"};

struct ReplayArm
{
  const char * name;
  Six old_total;
  Six new_total;
  Six saved_old_total;
  Six saved_new_total;
  StateValues old_state;
  Six expected_stress;
  StateValues expected_state;
  double expected_combined_damage;
};

State
unpackState(const StateValues & values)
{
  State state;
  std::copy_n(values.begin(), 6, state.backbone.plastic_strain.begin());
  state.backbone.tensile_equivalent_plastic_strain = values[6];
  state.backbone.compressive_equivalent_plastic_strain = values[7];
  std::copy_n(values.begin() + 8, 6, state.viscous_plastic_strain.begin());
  state.viscous_tension_damage = values[14];
  state.viscous_compression_damage = values[15];
  return state;
}

StateValues
packState(const State & state)
{
  StateValues values;
  std::copy(state.backbone.plastic_strain.begin(), state.backbone.plastic_strain.end(), values.begin());
  values[6] = state.backbone.tensile_equivalent_plastic_strain;
  values[7] = state.backbone.compressive_equivalent_plastic_strain;
  std::copy(state.viscous_plastic_strain.begin(), state.viscous_plastic_strain.end(), values.begin() + 8);
  values[14] = state.viscous_tension_damage;
  values[15] = state.viscous_compression_damage;
  return values;
}

void
expectTensorNear(const Six & actual,
                 const Six & expected,
                 const double tolerance,
                 const std::string & context)
{
  for (std::size_t i = 0; i < actual.size(); ++i)
  {
    SCOPED_TRACE(context + "." + tensor_names[i]);
    EXPECT_NEAR(actual[i], expected[i], tolerance);
  }
}

void
expectStateNear(const StateValues & actual,
                const StateValues & expected,
                const double plastic_tolerance,
                const double damage_tolerance,
                const std::string & context)
{
  for (std::size_t i = 0; i < actual.size(); ++i)
  {
    SCOPED_TRACE(context + "." + state_names[i]);
    EXPECT_NEAR(actual[i], expected[i], i < 14 ? plastic_tolerance : damage_tolerance);
  }
}

void
printReplay(const ReplayArm & arm, const AbaqusCDPSubstepIntegrator::LinearizedResult & replay)
{
  const auto & result = replay.result;
  const auto & final = result.final_result;
  const auto state = packState(final.state);
  std::cout << std::setprecision(17) << "DS04G_PEAK_REPLAY arm=" << arm.name
            << " N=" << result.accepted_substeps << " cutbacks=" << result.cutback_count
            << " partitions=" << result.attempted_partitions
            << " branch=" << static_cast<unsigned int>(final.backbone.active_branch)
            << " plastic=" << final.backbone.plastic << '\n';
  for (std::size_t i = 0; i < 6; ++i)
    std::cout << "DS04G_PEAK_REPLAY arm=" << arm.name << " stress_" << tensor_names[i]
              << " actual=" << final.cauchy_stress[i] << " saved=" << arm.expected_stress[i]
              << " delta=" << final.cauchy_stress[i] - arm.expected_stress[i] << '\n';
  for (std::size_t i = 0; i < state.size(); ++i)
    std::cout << "DS04G_PEAK_REPLAY arm=" << arm.name << " " << state_names[i]
              << " actual=" << state[i] << " saved=" << arm.expected_state[i]
              << " delta=" << state[i] - arm.expected_state[i] << '\n';
  std::cout << "DS04G_PEAK_REPLAY arm=" << arm.name
            << " combined actual=" << final.damage.damage
            << " saved=" << arm.expected_combined_damage
            << " delta=" << final.damage.damage - arm.expected_combined_damage << '\n';
}
}

// Values come from native-replay-input.json and frames 0052/0053 of the two locked 4a9 probes.
TEST(DS04GPeakReloadReplay, NativeSingleStepMatchesBothSavedFrames)
{
  const std::array<ReplayArm, 2> arms = {
      ReplayArm{"hold",
                {0.00019998000263977267, -5.0001250040585375e-05, -2.676624113323589e-19, 0.0, 0.0, 0.0},
                {0.0002000050026394927, -5.0001250040585375e-05, -2.676624113323589e-19, 0.0, 0.0, 0.0},
                {0.00019998000263977, -5.0001250040585e-05, 0.0, 0.0, 0.0, 0.0},
                {0.00020000500263949, -5.0001250040585e-05, 0.0, 0.0, 0.0, 0.0},
                {8.6114121775811e-05, -3.3865893666938e-05, -9.8680631026728e-06, 0.0, 0.0, 0.0,
                 8.6114121783318e-05, 0.0,
                 8.611412177449e-05, -3.3865893666418e-05, -9.8680631025213e-06, 0.0, 0.0, 0.0,
                 0.74641447500888, 0.0},
                {942716.118716, 124266.79077094, 287945.15020714, 0.0, 0.0, 0.0},
                {8.6127591249359e-05, -3.3871192034707e-05, -9.8696081660256e-06, 0.0, 0.0, 0.0,
                 8.6127591256867e-05, 0.0,
                 8.6120856511925e-05, -3.3868542850563e-05, -9.8688356342734e-06, 0.0, 0.0, 0.0,
                 0.74643983773763, 0.0},
                0.74643983773763},
      ReplayArm{"rotate",
                {0.00019998000264127212, -5.000125003953067e-05, -2.676624113323589e-19,
                 -2.3121894920629528e-17, 0.0, 0.0},
                {0.0002000050026409089, -5.0001250039447465e-05, -2.676624113323589e-19,
                 -1.4823005609374558e-17, 0.0, 0.0},
                {0.00019998000264126784, -5.000125003952986e-05, 0.0,
                 -2.4479675692122907e-17, 0.0, 0.0},
                {0.00020000500264090465, -5.0001250039446645e-05, 0.0,
                 -1.61808703913989e-17, 0.0, 0.0},
                {8.6114121775811e-05, -3.3865893666938e-05, -9.8680631026728e-06, 0.0, 0.0, 0.0,
                 8.6114121783318e-05, 0.0,
                 8.611412177449e-05, -3.3865893666418e-05, -9.8680631025213e-06, 0.0, 0.0, 0.0,
                 0.74641447500888, 0.0},
                {942716.11816002, 124266.8068056, 287945.15347453, 1.5735184196526, 0.0, 0.0},
                {8.6127589827182e-05, -3.3871191475164e-05, -9.8696080030156e-06,
                 3.6078502757122e-14, 0.0, 0.0,
                 8.6127589834689e-05, 0.0,
                 8.6120855800836e-05, -3.3868542570791e-05, -9.8688355527684e-06,
                 1.8039251378561e-14, 0.0, 0.0,
                 0.7464398350597, 0.0},
                0.7464398350597}};

  CDPMaterialTable table(data_dir + "compression_hardening.csv",
                         data_dir + "compression_damage.csv",
                         data_dir + "tension_stiffening.csv",
                         data_dir + "tension_damage.csv",
                         youngs_modulus);
  AbaqusCDPLocalIntegrator local(table, {youngs_modulus, 0.2, 30.0, 0.1, 1.16, 0.6667});
  AbaqusCDPStateIntegrator state_integrator(local, {0.0, 1.0, 0.0005});
  AbaqusCDPSubstepIntegrator integrator(state_integrator, {256, 2.5e-5, 1.0e-8});

  const auto hold = integrator.integrateLinearized(
      arms[0].old_total, arms[0].new_total, 0.0005, unpackState(arms[0].old_state));
  const auto rotate = integrator.integrateLinearized(
      arms[1].old_total, arms[1].new_total, 0.0005, unpackState(arms[1].old_state));
  const AbaqusCDPSubstepIntegrator::LinearizedResult * replays[] = {&hold, &rotate};

  for (std::size_t arm_index = 0; arm_index < arms.size(); ++arm_index)
  {
    const auto & arm = arms[arm_index];
    const auto & replay = *replays[arm_index];
    const auto & result = replay.result;
    const auto & final = result.final_result;
    printReplay(arm, replay);

    expectTensorNear(arm.old_total, arm.saved_old_total, 2.0e-15, std::string(arm.name) + ".old_total");
    expectTensorNear(arm.new_total, arm.saved_new_total, 2.0e-15, std::string(arm.name) + ".new_total");
    SCOPED_TRACE(std::string(arm.name) + ".counters");
    EXPECT_EQ(result.accepted_substeps, 1u);
    EXPECT_EQ(result.cutback_count, 0u);
    EXPECT_EQ(result.attempted_partitions, 1u);
    EXPECT_EQ(static_cast<unsigned int>(final.backbone.active_branch), 1u);
    EXPECT_TRUE(final.backbone.plastic);

    expectTensorNear(final.cauchy_stress, arm.expected_stress, 0.01,
                     std::string(arm.name) + ".saved_frame_0053.stress");
    expectStateNear(packState(final.state), arm.expected_state, 1.0e-12, 1.0e-10,
                    std::string(arm.name) + ".saved_frame_0053.state");
    EXPECT_NEAR(final.damage.damage, arm.expected_combined_damage, 1.0e-10)
        << arm.name << ".saved_frame_0053.combined";
  }

  expectTensorNear(hold.result.final_result.cauchy_stress,
                   rotate.result.final_result.cauchy_stress, 0.01, "between_arms.stress");
  expectStateNear(packState(hold.result.final_result.state),
                  packState(rotate.result.final_result.state), 1.0e-13, 1.0e-9, "between_arms.state");
  EXPECT_NEAR(hold.result.final_result.damage.damage,
              rotate.result.final_result.damage.damage, 1.0e-9)
      << "between_arms.combined";
}
