#include "CDPDiagnostics.h"
#include "AbaqusCDPStateIntegrator.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace
{
using SymmetricTensor = AbaqusCDPStateIntegrator::SymmetricTensor;

[[noreturn]] void
stateIntegrationError(const std::string & message)
{
  throw std::runtime_error("AbaqusCDPStateIntegrator: " + message);
}

bool
finiteStateTensor(const SymmetricTensor & tensor)
{
  return std::all_of(
      tensor.begin(), tensor.end(), [](const double value) { return std::isfinite(value); });
}

SymmetricTensor
subtractStateTensor(const SymmetricTensor & left, const SymmetricTensor & right)
{
  SymmetricTensor result;
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = left[i] - right[i];
  return result;
}

SymmetricTensor
scaleStateTensor(const SymmetricTensor & tensor, const double factor)
{
  SymmetricTensor result;
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = factor * tensor[i];
  return result;
}

double
stateTensorNorm(const SymmetricTensor & tensor)
{
  double sum = 0.0;
  for (const double value : tensor)
    sum += value * value;
  return std::sqrt(sum);
}

void
validateDamage(const double damage, const std::string & name)
{
  if (!std::isfinite(damage) || damage < 0.0 || damage >= 1.0)
    stateIntegrationError(name + " must be finite and in [0, 1)");
}
}

AbaqusCDPStateIntegrator::AbaqusCDPStateIntegrator(
    const AbaqusCDPLocalIntegrator & backbone_integrator,
    Parameters parameters)
  : _backbone_integrator(backbone_integrator), _parameters(parameters)
{
  if (!std::isfinite(_parameters.tension_recovery) ||
      !std::isfinite(_parameters.compression_recovery) ||
      _parameters.tension_recovery < 0.0 || _parameters.tension_recovery > 1.0 ||
      _parameters.compression_recovery < 0.0 || _parameters.compression_recovery > 1.0)
    stateIntegrationError("recovery factors must be finite and in [0, 1]");
  if (!std::isfinite(_parameters.relaxation_time) || _parameters.relaxation_time < 0.0)
    stateIntegrationError("relaxation time must be finite and nonnegative");
  if (!std::isfinite(_parameters.state_tolerance) || _parameters.state_tolerance < 0.0)
    stateIntegrationError("state tolerance must be finite and nonnegative");
}

AbaqusCDPStateIntegrator::Result
AbaqusCDPStateIntegrator::integrate(const SymmetricTensor & total_strain,
                                    const double time_step,
                                    const State & old_state) const
{
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::STATE);
  if (!finiteStateTensor(total_strain) || !finiteStateTensor(old_state.viscous_plastic_strain))
    stateIntegrationError("total strain or old viscous plastic strain contains a non-finite value");
  if (!std::isfinite(time_step) || time_step < 0.0)
    stateIntegrationError("time step must be finite and nonnegative");
  validateDamage(old_state.viscous_tension_damage, "old viscous tension damage");
  validateDamage(old_state.viscous_compression_damage, "old viscous compression damage");
  validateDamage(old_state.viscous_combined_damage, "old viscous combined damage");

  // The old state remains immutable. Any failure below leaves the caller's
  // checkpoint untouched and makes a cutback retry deterministic.
  const auto backbone = _backbone_integrator.integrate(total_strain, old_state.backbone);
  return assembleResult(total_strain, time_step, old_state, backbone);
}

AbaqusCDPStateIntegrator::Result
AbaqusCDPStateIntegrator::assembleResult(
    const SymmetricTensor & total_strain,
    const double time_step,
    const State & old_state,
    const AbaqusCDPLocalIntegrator::Result & backbone) const
{
  validateDamage(backbone.backbone_tension_damage, "backbone tension damage");
  validateDamage(backbone.backbone_compression_damage, "backbone compression damage");
  if (backbone.backbone_tension_damage + _parameters.state_tolerance <
      old_state.viscous_tension_damage)
    stateIntegrationError("backbone tension damage is below old viscous tension damage");
  if (backbone.backbone_compression_damage + _parameters.state_tolerance <
      old_state.viscous_compression_damage)
    stateIntegrationError("backbone compression damage is below old viscous compression damage");

  State new_state;
  new_state.backbone = backbone.state;
  new_state.viscous_plastic_strain = AbaqusCDPFormula::duvautLionsUpdate(
      old_state.viscous_plastic_strain,
      backbone.state.plastic_strain,
      time_step,
      _parameters.relaxation_time);
  new_state.viscous_tension_damage = AbaqusCDPFormula::duvautLionsUpdate(
      old_state.viscous_tension_damage,
      backbone.backbone_tension_damage,
      time_step,
      _parameters.relaxation_time);
  new_state.viscous_compression_damage = AbaqusCDPFormula::duvautLionsUpdate(
      old_state.viscous_compression_damage,
      backbone.backbone_compression_damage,
      time_step,
      _parameters.relaxation_time);

  validateDamage(new_state.viscous_tension_damage, "new viscous tension damage");
  validateDamage(new_state.viscous_compression_damage, "new viscous compression damage");
  if (new_state.viscous_tension_damage + _parameters.state_tolerance <
      old_state.viscous_tension_damage ||
      new_state.viscous_compression_damage + _parameters.state_tolerance <
          old_state.viscous_compression_damage)
    stateIntegrationError("viscous damage must be irreversible");

  const double backbone_tension_weight =
      AbaqusCDPFormula::stressInvariants(backbone.effective_stress).tension_weight;
  const auto backbone_damage = AbaqusCDPFormula::combineDamage(
      backbone.backbone_compression_damage,
      backbone.backbone_tension_damage,
      _parameters.tension_recovery,
      _parameters.compression_recovery,
      backbone_tension_weight);
  new_state.viscous_combined_damage =
      _parameters.use_scalar_damage_viscosity
          ? AbaqusCDPFormula::duvautLionsUpdate(old_state.viscous_combined_damage,
                                                backbone_damage.damage,
                                                time_step,
                                                _parameters.relaxation_time)
          : 0.0;
  validateDamage(new_state.viscous_combined_damage, "new viscous combined damage");

  const auto viscous_effective_stress = _backbone_integrator.elasticStress(
      subtractStateTensor(total_strain, new_state.viscous_plastic_strain));
  const double tension_weight =
      AbaqusCDPFormula::stressInvariants(viscous_effective_stress).tension_weight;
  auto damage = AbaqusCDPFormula::combineDamage(new_state.viscous_compression_damage,
                                                new_state.viscous_tension_damage,
                                                _parameters.tension_recovery,
                                                _parameters.compression_recovery,
                                                tension_weight);
  if (_parameters.use_scalar_damage_viscosity)
  {
    damage.damage = new_state.viscous_combined_damage;
    damage.stiffness_factor = 1.0 - new_state.viscous_combined_damage;
  }
  else
    new_state.viscous_combined_damage = damage.damage;
  const auto cauchy_stress = scaleStateTensor(viscous_effective_stress, damage.stiffness_factor);

  const auto plastic_lag =
      subtractStateTensor(backbone.state.plastic_strain, new_state.viscous_plastic_strain);
  const bool rate_independent = _parameters.relaxation_time == 0.0;
  const double dt_over_relaxation_time =
      rate_independent ? 0.0 : time_step / _parameters.relaxation_time;

  return {backbone,
          viscous_effective_stress,
          cauchy_stress,
          new_state,
          damage,
          rate_independent,
          dt_over_relaxation_time,
          stateTensorNorm(plastic_lag),
          backbone.backbone_tension_damage - new_state.viscous_tension_damage,
          backbone.backbone_compression_damage - new_state.viscous_compression_damage,
          backbone_damage.damage,
          backbone_damage.damage - new_state.viscous_combined_damage};
}

AbaqusCDPStateIntegrator::LinearizedResult
AbaqusCDPStateIntegrator::integrateLinearized(const SymmetricTensor & total_strain,
                                              const double time_step,
                                              const State & old_state) const
{
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::STATE_LINEARIZED);
  if (!finiteStateTensor(total_strain) || !finiteStateTensor(old_state.viscous_plastic_strain))
    stateIntegrationError("total strain or old viscous plastic strain contains a non-finite value");
  if (!std::isfinite(time_step) || time_step < 0.0)
    stateIntegrationError("time step must be finite and nonnegative");
  validateDamage(old_state.viscous_tension_damage, "old viscous tension damage");
  validateDamage(old_state.viscous_compression_damage, "old viscous compression damage");
  validateDamage(old_state.viscous_combined_damage, "old viscous combined damage");

  const auto backbone =
      _backbone_integrator.integrateLinearized(total_strain, old_state.backbone);
  LinearizedResult linearized{
      assembleResult(total_strain, time_step, old_state, backbone.result), {}};

  const double relaxation_factor =
      _parameters.relaxation_time == 0.0
          ? 1.0
          : (time_step == 0.0 ? 0.0
                              : time_step / (_parameters.relaxation_time + time_step));
  const auto tension = _backbone_integrator.materialResponse(
      CDPMaterialTable::Branch::TENSION,
      linearized.result.state.backbone.tensile_equivalent_plastic_strain);
  const auto compression = _backbone_integrator.materialResponse(
      CDPMaterialTable::Branch::COMPRESSION,
      linearized.result.state.backbone.compressive_equivalent_plastic_strain);

  const auto weightGradient = [](const SymmetricTensor & stress) {
    SymmetricTensor gradient;
    double stress_scale = 1.0;
    for (const double value : stress)
      stress_scale = std::max(stress_scale, std::abs(value));
    const double weight_step = 1.0e-7 * stress_scale;
    for (std::size_t column = 0; column < 6; ++column)
    {
      auto plus = stress;
      auto minus = stress;
      plus[column] += weight_step;
      minus[column] -= weight_step;
      gradient[column] =
          (AbaqusCDPFormula::stressInvariants(plus).tension_weight -
           AbaqusCDPFormula::stressInvariants(minus).tension_weight) /
          (2.0 * weight_step);
    }
    return gradient;
  };
  const auto tension_weight_gradient = weightGradient(linearized.result.viscous_effective_stress);
  const auto backbone_weight_gradient = weightGradient(linearized.result.backbone.effective_stress);

  const auto backbone_damage = AbaqusCDPFormula::combineDamage(
      linearized.result.backbone.backbone_compression_damage,
      linearized.result.backbone.backbone_tension_damage,
      _parameters.tension_recovery,
      _parameters.compression_recovery,
      AbaqusCDPFormula::stressInvariants(linearized.result.backbone.effective_stress)
          .tension_weight);
  const double backbone_compression_factor =
      1.0 - backbone_damage.tensile_recovery_factor *
                linearized.result.backbone.backbone_compression_damage;
  const double backbone_tension_factor =
      1.0 - backbone_damage.compressive_recovery_factor *
                linearized.result.backbone.backbone_tension_damage;

  const double tension_recovery_factor = linearized.result.damage.tensile_recovery_factor;
  const double compression_recovery_factor =
      linearized.result.damage.compressive_recovery_factor;
  const double compression_factor =
      1.0 - tension_recovery_factor * linearized.result.state.viscous_compression_damage;
  const double tension_factor =
      1.0 - compression_recovery_factor * linearized.result.state.viscous_tension_damage;

  for (std::size_t input = 0; input < transition_size; ++input)
  {
    std::array<double, AbaqusCDPLocalIntegrator::transition_size> local_input = {};
    if (input < 6)
      local_input[input] = 1.0;
    else if (input - 6 < 8)
      local_input[6 + input - 6] = 1.0;

    std::array<double, AbaqusCDPLocalIntegrator::transition_size> local_output = {};
    for (std::size_t local_column = 0;
         local_column < AbaqusCDPLocalIntegrator::transition_size;
         ++local_column)
      if (local_input[local_column] != 0.0)
        for (std::size_t row = 0; row < AbaqusCDPLocalIntegrator::transition_size; ++row)
          local_output[row] +=
              backbone.derivative[local_column][row] * local_input[local_column];

    SymmetricTensor strain_derivative = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (input < 6)
      strain_derivative[input] = 1.0;
    SymmetricTensor old_viscous_plastic_derivative = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double old_tension_damage_derivative = 0.0;
    double old_compression_damage_derivative = 0.0;
    double old_combined_damage_derivative = 0.0;
    if (input >= 6)
    {
      const std::size_t old_state_index = input - 6;
      if (old_state_index >= 8 && old_state_index < 14)
        old_viscous_plastic_derivative[old_state_index - 8] = 1.0;
      else if (old_state_index == 14)
        old_tension_damage_derivative = 1.0;
      else if (old_state_index == 15)
        old_compression_damage_derivative = 1.0;
      else if (old_state_index == 16)
        old_combined_damage_derivative = 1.0;
    }

    SymmetricTensor new_viscous_plastic_derivative;
    for (std::size_t row = 0; row < 6; ++row)
      new_viscous_plastic_derivative[row] =
          (1.0 - relaxation_factor) * old_viscous_plastic_derivative[row] +
          relaxation_factor * local_output[6 + row];
    const double new_tension_damage_derivative =
        (1.0 - relaxation_factor) * old_tension_damage_derivative +
        relaxation_factor * tension.damage.right_derivative * local_output[12];
    const double new_compression_damage_derivative =
        (1.0 - relaxation_factor) * old_compression_damage_derivative +
        relaxation_factor * compression.damage.right_derivative * local_output[13];

    SymmetricTensor elastic_argument_derivative;
    for (std::size_t row = 0; row < 6; ++row)
      elastic_argument_derivative[row] =
          strain_derivative[row] - new_viscous_plastic_derivative[row];
    const auto effective_stress_derivative =
        _backbone_integrator.elasticStress(elastic_argument_derivative);
    double tension_weight_derivative = 0.0;
    double backbone_weight_derivative = 0.0;
    for (std::size_t row = 0; row < 6; ++row)
    {
      tension_weight_derivative +=
          tension_weight_gradient[row] * effective_stress_derivative[row];
      backbone_weight_derivative += backbone_weight_gradient[row] * local_output[row];
    }

    const double tension_recovery_derivative =
        -_parameters.tension_recovery * tension_weight_derivative;
    const double compression_recovery_derivative =
        _parameters.compression_recovery * tension_weight_derivative;
    const double compression_factor_derivative =
        -(tension_recovery_derivative *
              linearized.result.state.viscous_compression_damage +
          tension_recovery_factor * new_compression_damage_derivative);
    const double tension_factor_derivative =
        -(compression_recovery_derivative *
              linearized.result.state.viscous_tension_damage +
          compression_recovery_factor * new_tension_damage_derivative);
    const double branch_stiffness_derivative =
        compression_factor_derivative * tension_factor +
        compression_factor * tension_factor_derivative;

    const double backbone_tension_recovery_derivative =
        -_parameters.tension_recovery * backbone_weight_derivative;
    const double backbone_compression_recovery_derivative =
        _parameters.compression_recovery * backbone_weight_derivative;
    const double backbone_compression_factor_derivative =
        -(backbone_tension_recovery_derivative *
              linearized.result.backbone.backbone_compression_damage +
          backbone_damage.tensile_recovery_factor * compression.damage.right_derivative *
              local_output[13]);
    const double backbone_tension_factor_derivative =
        -(backbone_compression_recovery_derivative *
              linearized.result.backbone.backbone_tension_damage +
          backbone_damage.compressive_recovery_factor * tension.damage.right_derivative *
              local_output[12]);
    const double backbone_damage_derivative =
        -(backbone_compression_factor_derivative * backbone_tension_factor +
          backbone_compression_factor * backbone_tension_factor_derivative);
    const double new_combined_damage_derivative =
        _parameters.use_scalar_damage_viscosity
            ? (1.0 - relaxation_factor) * old_combined_damage_derivative +
                  relaxation_factor * backbone_damage_derivative
            : -branch_stiffness_derivative;
    const double stiffness_derivative =
        _parameters.use_scalar_damage_viscosity ? -new_combined_damage_derivative
                                                : branch_stiffness_derivative;

    for (std::size_t row = 0; row < 6; ++row)
      linearized.derivative[input][row] =
          linearized.result.damage.stiffness_factor * effective_stress_derivative[row] +
          stiffness_derivative * linearized.result.viscous_effective_stress[row];
    for (std::size_t row = 0; row < 6; ++row)
    {
      linearized.derivative[input][6 + row] = local_output[6 + row];
      linearized.derivative[input][6 + 8 + row] =
          new_viscous_plastic_derivative[row];
    }
    linearized.derivative[input][6 + 6] = local_output[12];
    linearized.derivative[input][6 + 7] = local_output[13];
    linearized.derivative[input][6 + 14] = new_tension_damage_derivative;
    linearized.derivative[input][6 + 15] = new_compression_damage_derivative;
    linearized.derivative[input][6 + 16] = new_combined_damage_derivative;
  }
  return linearized;
}
