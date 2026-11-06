#include "AbaqusCDPSubstepIntegrator.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

namespace
{
struct AcceptedPathFrame
{
  char stage;
  int increment;
  double time_old;
  double time_new;
  double dt;
  Six D;
};

struct MidProxy
{
  std::array<double, 4> stress;
};

struct ReplayArmSummary
{
  int first_plastic_frame = 0;
  std::vector<unsigned int> branches;
  double final_time = 0.0;
};

void
requireReplay(bool condition, const std::string & message)
{
  if (!condition)
    throw std::runtime_error(message);
}

std::vector<AcceptedPathFrame>
readAcceptedPath()
{
  std::ifstream input(data_dir + "959_v2_3_accepted_path.txt");
  requireReplay(input.is_open(), "cannot open 959 accepted path");
  std::string identity;
  unsigned int version;
  unsigned int count;
  input >> identity >> version >> count;
  requireReplay(identity == "DS04G_959_V2_3_PATH" && version == 1 && count == 187,
                "unexpected 959 path identity/version/count");
  std::vector<AcceptedPathFrame> frames;
  for (unsigned int index = 0; index < count; ++index)
  {
    AcceptedPathFrame frame;
    input >> frame.stage >> frame.increment >> frame.time_old >> frame.time_new >> frame.dt;
    for (auto & value : frame.D)
      input >> value;
    requireReplay(input.good(), "truncated 959 path at frame " + std::to_string(index + 1));
    requireReplay(frame.stage == (index < 24 ? 'G' : 'E'), "stage order differs from G24/E163");
    requireReplay(frame.increment == static_cast<int>(index < 24 ? index + 1 : index - 23),
                  "accepted increment sequence mismatch");
    frames.push_back(frame);
  }
  std::string extra;
  requireReplay(!(input >> extra), "extra content after 187 accepted increments");
  return frames;
}

std::vector<std::string>
csvColumns(const std::string & line)
{
  std::vector<std::string> columns;
  std::stringstream stream(line);
  std::string column;
  while (std::getline(stream, column, ','))
    columns.push_back(column);
  return columns;
}

std::map<std::string, MidProxy>
readMidProxy()
{
  std::ifstream input(data_dir + "959_v2_3_mid_proxy.csv");
  requireReplay(input.is_open(), "cannot open 959 M_mid proxy CSV");
  std::string line;
  requireReplay(static_cast<bool>(std::getline(input, line)), "empty 959 reference CSV");
  const auto header = csvColumns(line);
  std::map<std::string, std::size_t> column;
  for (std::size_t i = 0; i < header.size(); ++i)
    column[header[i]] = i;
  const std::array<const char *, 4> names = {"S11", "S22", "S33", "S12"};
  for (const auto key : {"step", "increment"})
    requireReplay(column.count(key), "reference missing key column " + std::string(key));
  for (const auto * name : names)
    requireReplay(column.count(std::string("midpoint_cayley_S_Pa_") + name),
                  "reference missing M_mid component " + std::string(name));

  std::map<std::string, MidProxy> result;
  while (std::getline(input, line))
  {
    if (line.empty())
      continue;
    const auto values = csvColumns(line);
    requireReplay(values.size() == header.size(), "reference CSV row width mismatch");
    const auto & step = values[column["step"]];
    const char stage = step == "Gravity" ? 'G' : step == "Earthquake" ? 'E' : '?';
    requireReplay(stage != '?', "unknown reference step " + step);
    MidProxy proxy;
    for (std::size_t i = 0; i < names.size(); ++i)
    {
      proxy.stress[i] = std::stod(values[column[std::string("midpoint_cayley_S_Pa_") + names[i]]]);
      requireReplay(std::isfinite(proxy.stress[i]), "nonfinite M_mid proxy");
    }
    const auto key = std::string(1, stage) + ":" + values[column["increment"]];
    requireReplay(!result.count(key), "duplicate reference key " + key);
    result[key] = proxy;
  }
  return result;
}

double
float32Ulp(double value)
{
  const float single = static_cast<float>(value);
  const float following = std::nextafter(single, std::numeric_limits<float>::infinity());
  return std::abs(static_cast<double>(following) - static_cast<double>(single));
}

std::string
jsonNumber(double value)
{
  if (!std::isfinite(value))
    return "null";
  std::ostringstream output;
  output << std::setprecision(17) << value;
  return output.str();
}

template <std::size_t N>
std::string
jsonArray(const std::array<double, N> & values)
{
  std::string result = "[";
  for (std::size_t i = 0; i < values.size(); ++i)
  {
    if (i)
      result += ',';
    result += jsonNumber(values[i]);
  }
  return result + ']';
}

std::string
jsonEscape(const std::string & value)
{
  std::string result;
  for (const unsigned char c : value)
  {
    if (c == '"' || c == '\\')
      result += '\\', result += static_cast<char>(c);
    else if (c == '\n')
      result += "\\n";
    else if (c == '\r')
      result += "\\r";
    else if (c == '\t')
      result += "\\t";
    else if (c >= 0x20)
      result += static_cast<char>(c);
  }
  return result;
}

ReplayArmSummary
run959Arm(const std::string & arm,
          int dt_direction,
          const std::vector<AcceptedPathFrame> & path,
          const std::map<std::string, MidProxy> & proxies,
          const AbaqusCDPSubstepIntegrator & integrator)
{
  State old_state{};
  Six old_total{};
  Six elastic_stress{};
  double current_time = 0.0;
  std::size_t used_proxies = 0;
  bool elastic_gate = true;
  ReplayArmSummary summary;
  const double lambda = youngs_modulus * 0.2 / (1.2 * 0.6);
  const double mu = youngs_modulus / 2.4;

  for (std::size_t index = 0; index < path.size(); ++index)
  {
    const auto & frame = path[index];
    const int frame_number = static_cast<int>(index + 1);
    const bool eq_restart = frame.stage == 'E' && frame.increment == 1;
    requireReplay(frame.time_old == (eq_restart ? 0.0 : current_time),
                  arm + ": noncontiguous accepted time outside the Gravity/EQ handoff");
    requireReplay(frame.time_new > frame.time_old && frame.dt > 0,
                  arm + ": nonpositive source dt");
    requireReplay(frame.dt == frame.time_new - frame.time_old,
                  arm + ": source dt differs from native frame time difference");
    for (const double value : frame.D)
      requireReplay(std::isfinite(value), arm + ": nonfinite D");
    requireReplay(frame.D[2] == 0 && frame.D[4] == 0 && frame.D[5] == 0,
                  arm + ": path is not the frozen planar D representation");

    const double clock_delta = float32Ulp(frame.time_new) + float32Ulp(frame.time_old) + 2e-9;
    const double dt = dt_direction == 0
                          ? frame.dt
                          : std::max(0.0, frame.dt + dt_direction * clock_delta);
    requireReplay(dt > 0, arm + ": nonpositive sensitivity dt");

    Six new_total;
    Six recovered;
    for (std::size_t component = 0; component < 6; ++component)
    {
      new_total[component] = old_total[component] + frame.D[component];
      recovered[component] = new_total[component] - old_total[component];
      const double scale = std::max({1.0,
                                     std::abs(old_total[component]),
                                     std::abs(new_total[component]),
                                     std::abs(frame.D[component])});
      const double budget = 8 * std::numeric_limits<double>::epsilon() * scale;
      requireReplay(std::abs(recovered[component] - frame.D[component]) <= budget,
                    arm + ": old/new total-D rounding budget exceeded");
    }

    Six elastic_candidate = elastic_stress;
    const double trace = frame.D[0] + frame.D[1] + frame.D[2];
    for (std::size_t component = 0; component < 6; ++component)
    {
      elastic_candidate[component] += 2 * mu * frame.D[component];
      if (component < 3)
        elastic_candidate[component] += lambda * trace;
    }

    const auto old_state_before = old_state;
    AbaqusCDPSubstepIntegrator::LinearizedResult replay;
    try
    {
      replay = integrator.integrateLinearized(old_total, new_total, dt, old_state);
    }
    catch (const std::exception & error)
    {
      std::cout << "{\"record_type\":\"native_exception\",\"arm\":\"" << arm
                << "\",\"frame\":" << frame_number << ",\"time_s\":" << jsonNumber(frame.time_new)
                << ",\"dt_s\":" << jsonNumber(dt) << ",\"error\":\"" << jsonEscape(error.what())
                << "\"}\n";
      throw;
    }
    catch (...)
    {
      std::cout << "{\"record_type\":\"native_exception\",\"arm\":\"" << arm
                << "\",\"frame\":" << frame_number << ",\"error\":\"non-standard exception\"}\n";
      throw;
    }

    const auto & result = replay.result;
    const auto & final = result.final_result;
    const auto state16 = packState(final.state);
    std::string proxy_stress = "null";
    std::string proxy_delta = "null";
    const auto proxy = proxies.find(std::string(1, frame.stage) + ":" +
                                    std::to_string(frame.increment));
    if (proxy != proxies.end())
    {
      const std::array<double, 4> native_four = {final.cauchy_stress[0],
                                                 final.cauchy_stress[1],
                                                 final.cauchy_stress[2],
                                                 final.cauchy_stress[3]};
      std::array<double, 4> delta;
      for (std::size_t i = 0; i < 4; ++i)
        delta[i] = native_four[i] - proxy->second.stress[i];
      proxy_stress = jsonArray(proxy->second.stress);
      proxy_delta = jsonArray(delta);
      ++used_proxies;
    }
    Six elastic_delta;
    for (std::size_t i = 0; i < 6; ++i)
      elastic_delta[i] = final.cauchy_stress[i] - elastic_candidate[i];
    const bool old_unchanged = packState(old_state) == packState(old_state_before);

    std::cout << "{\"record_type\":\"frame\",\"arm\":\"" << arm << "\",\"frame\":" << frame_number
              << ",\"stage\":\"" << frame.stage << "\",\"increment\":" << frame.increment
              << ",\"time_s\":" << jsonNumber(frame.time_new) << ",\"dt_s\":" << jsonNumber(dt)
              << ",\"clock_delta_s\":" << jsonNumber(clock_delta) << ",\"D\":" << jsonArray(frame.D)
              << ",\"old_total\":" << jsonArray(old_total) << ",\"new_total\":" << jsonArray(new_total)
              << ",\"cauchy_stress\":" << jsonArray(final.cauchy_stress)
              << ",\"state16\":" << jsonArray(state16)
              << ",\"combined_damage\":" << jsonNumber(final.damage.damage)
              << ",\"stiffness_factor\":" << jsonNumber(final.damage.stiffness_factor)
              << ",\"branch\":" << static_cast<unsigned int>(final.backbone.active_branch)
              << ",\"plastic\":" << (final.backbone.plastic ? 1 : 0)
              << ",\"N\":" << result.accepted_substeps
              << ",\"cutbacks\":" << result.cutback_count
              << ",\"partitions\":" << result.attempted_partitions
              << ",\"final_local_residual\":" << jsonNumber(final.backbone.residual_norm)
              << ",\"old_state_unchanged\":" << (old_unchanged ? 1 : 0)
              << ",\"elastic_reference_stress\":" << jsonArray(elastic_candidate)
              << ",\"elastic_delta\":" << jsonArray(elastic_delta)
              << ",\"elastic_gate_applied\":" << (elastic_gate ? 1 : 0)
              << ",\"mid_proxy_stress\":" << proxy_stress
              << ",\"mid_proxy_delta\":" << proxy_delta << "}\n";

    requireReplay(old_unchanged, arm + ": integrator modified old_state");
    requireReplay(result.accepted_substeps > 0 && result.accepted_substeps <= 256,
                  arm + ": accepted substeps outside 1..256");
    requireReplay(result.attempted_partitions > 0, arm + ": missing native partition");
    for (const double value : final.cauchy_stress)
      requireReplay(std::isfinite(value), arm + ": nonfinite Cauchy stress");
    for (const double value : state16)
      requireReplay(std::isfinite(value), arm + ": nonfinite state16");
    requireReplay(std::isfinite(final.damage.damage) && final.damage.damage >= 0 &&
                      final.damage.damage <= 1,
                  arm + ": combined damage outside [0,1]");
    requireReplay(std::isfinite(final.damage.stiffness_factor) &&
                      final.damage.stiffness_factor >= 0 && final.damage.stiffness_factor <= 1,
                  arm + ": stiffness factor outside [0,1]");
    requireReplay(std::isfinite(final.backbone.residual_norm) && final.backbone.residual_norm >= 0,
                  arm + ": invalid final local residual");
    requireReplay(state16[6] >= packState(old_state_before)[6] &&
                      state16[7] >= packState(old_state_before)[7] &&
                      state16[14] >= packState(old_state_before)[14] &&
                      state16[15] >= packState(old_state_before)[15],
                  arm + ": kappa or branch damage decreased");

    if (elastic_gate)
    {
      for (std::size_t component = 0; component < 6; ++component)
      {
        const double scale = std::max(std::abs(final.cauchy_stress[component]),
                                      std::abs(elastic_candidate[component]));
        const double budget = 1e-5 + 1e-12 * scale;
        requireReplay(std::abs(elastic_delta[component]) <= budget,
                      arm + ": independent unrotated elastic recurrence mismatch at component " +
                          std::to_string(component));
      }
    }
    if (final.backbone.plastic)
    {
      elastic_gate = false;
      if (!summary.first_plastic_frame)
        summary.first_plastic_frame = frame_number;
    }
    summary.branches.push_back(static_cast<unsigned int>(final.backbone.active_branch));
    old_state = final.state;
    old_total = new_total;
    elastic_stress = elastic_candidate;
    current_time = frame.time_new;
  }
  requireReplay(used_proxies == 118, arm + ": M_mid elastic-window coverage mismatch");
  summary.final_time = current_time;
  requireReplay(current_time == 2.7999999523162842,
                arm + ": replay did not reach the frozen EQ2.80 frame");
  std::cout << "{\"record_type\":\"arm_summary\",\"arm\":\"" << arm
            << "\",\"frames\":187,\"first_plastic_frame\":" << summary.first_plastic_frame
            << ",\"final_time_s\":" << jsonNumber(summary.final_time) << "}\n";
  return summary;
}
}

TEST(DS04GPeakReloadReplay, Native959AcceptedPath)
{
  try
  {
    const auto path = readAcceptedPath();
    const auto proxies = readMidProxy();
    requireReplay(path.size() == 187, "959 path must contain 187 increments");

    CDPMaterialTable table(data_dir + "compression_hardening.csv",
                           data_dir + "compression_damage.csv",
                           data_dir + "tension_stiffening.csv",
                           data_dir + "tension_damage.csv",
                           youngs_modulus);
    AbaqusCDPLocalIntegrator local(table, {youngs_modulus, 0.2, 30.0, 0.1, 1.16, 0.6667});
    AbaqusCDPStateIntegrator state_integrator(local, {0.0, 1.0, 0.0005});
    AbaqusCDPSubstepIntegrator integrator(state_integrator, {256, 2.5e-5, 1.0e-8});

    const auto native = run959Arm("native", 0, path, proxies, integrator);
    const auto plus = run959Arm("dt_plus", 1, path, proxies, integrator);
    const auto minus = run959Arm("dt_minus", -1, path, proxies, integrator);
    const bool clock_sensitive = native.first_plastic_frame != plus.first_plastic_frame ||
                                 native.first_plastic_frame != minus.first_plastic_frame ||
                                 native.branches != plus.branches || native.branches != minus.branches;
    std::cout << "{\"record_type\":\"clock_sensitivity\",\"arms\":3,"
              << "\"first_plastic_frames\":[" << native.first_plastic_frame << ','
              << plus.first_plastic_frame << ',' << minus.first_plastic_frame << "],"
              << "\"branch_sequences_equal\":"
              << (native.branches == plus.branches && native.branches == minus.branches ? 1 : 0)
              << ",\"clock_sensitive\":" << (clock_sensitive ? 1 : 0)
              << ",\"strict_clock_closed\":0}\n";
  }
  catch (const std::exception & error)
  {
    std::cout << "{\"record_type\":\"hard_failure\",\"error\":\"" << jsonEscape(error.what())
              << "\"}\n";
    FAIL() << "959 native accepted-path replay failed: " << error.what();
  }
  catch (...)
  {
    std::cout << "{\"record_type\":\"hard_failure\",\"error\":\"non-standard exception\"}\n";
    FAIL() << "959 native accepted-path replay failed with a non-standard exception";
  }
}
