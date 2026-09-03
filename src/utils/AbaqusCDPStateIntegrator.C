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
  validateAssemblyInputs(total_strain, time_step, old_state);

  // The old state remains immutable. Any failure below leaves the caller's
  // checkpoint untouched and makes a cutback retry deterministic.
  const auto backbone = integrateBackbone(total_strain, old_state.backbone);
  return assembleResult(total_strain, time_step, old_state, backbone);
}

void
AbaqusCDPStateIntegrator::validateAssemblyInputs(const SymmetricTensor & total_strain,
                                                 const double time_step,
                                                 const State & old_state) const
{
  if (!finiteStateTensor(total_strain) || !finiteStateTensor(old_state.viscous_plastic_strain))
    stateIntegrationError("total strain or old viscous plastic strain contains a non-finite value");
  if (!std::isfinite(time_step) || time_step < 0.0)
    stateIntegrationError("time step must be finite and nonnegative");
  validateDamage(old_state.viscous_tension_damage, "old viscous tension damage");
  validateDamage(old_state.viscous_compression_damage, "old viscous compression damage");
}

AbaqusCDPLocalIntegrator::Result
AbaqusCDPStateIntegrator::integrateBackbone(
    const SymmetricTensor & total_strain,
    const AbaqusCDPLocalIntegrator::State & old_state) const
{
  return _backbone_integrator.integrate(total_strain, old_state);
}

AbaqusCDPLocalIntegrator::LinearizedResult
AbaqusCDPStateIntegrator::integrateBackboneLinearized(
    const SymmetricTensor & total_strain,
    const AbaqusCDPLocalIntegrator::State & old_state) const
{
  return _backbone_integrator.integrateLinearized(total_strain, old_state);
}

AbaqusCDPLocalIntegrator::Result
AbaqusCDPStateIntegrator::aggregateBackboneHistory(
    const AbaqusCDPLocalIntegrator::State & old_backbone_state,
    const AbaqusCDPLocalIntegrator::Result & substepped_backbone) const
{
  const double tension_weight =
      AbaqusCDPFormula::stressInvariants(substepped_backbone.effective_stress).tension_weight;
  return aggregateBackboneHistory(
      old_backbone_state, substepped_backbone, tension_weight, 1.0 - tension_weight);
}

AbaqusCDPLocalIntegrator::Result
AbaqusCDPStateIntegrator::aggregateBackboneHistory(
    const AbaqusCDPLocalIntegrator::State & old_backbone_state,
    const AbaqusCDPLocalIntegrator::Result & substepped_backbone,
    const double tensile_path_weight,
    const double compressive_path_weight) const
{
  if (!std::isfinite(tensile_path_weight) || tensile_path_weight < 0.0 ||
      tensile_path_weight > 1.0 || !std::isfinite(compressive_path_weight) ||
      compressive_path_weight < 0.0 || compressive_path_weight > 1.0)
    stateIntegrationError("accepted-step history path weights must be finite and in [0, 1]");

  auto aggregated = substepped_backbone;
  const auto plastic_increment = subtractStateTensor(
      substepped_backbone.state.plastic_strain, old_backbone_state.plastic_strain);
  const auto plastic_principal =
      AbaqusCDPFormula::stressInvariants(plastic_increment).principal_stress;
  const double tensile_increment =
      tensile_path_weight * std::max(0.0, plastic_principal.back());
  const double compressive_increment =
      compressive_path_weight * std::max(0.0, -plastic_principal.front());

  aggregated.state.tensile_equivalent_plastic_strain =
      old_backbone_state.tensile_equivalent_plastic_strain + tensile_increment;
  aggregated.state.compressive_equivalent_plastic_strain =
      old_backbone_state.compressive_equivalent_plastic_strain + compressive_increment;
  aggregated.backbone_tension_damage =
      _backbone_integrator
          .materialResponse(CDPMaterialTable::Branch::TENSION,
                            aggregated.state.tensile_equivalent_plastic_strain)
          .damage.value;
  aggregated.backbone_compression_damage =
      _backbone_integrator
          .materialResponse(CDPMaterialTable::Branch::COMPRESSION,
                            aggregated.state.compressive_equivalent_plastic_strain)
          .damage.value;

  return aggregated;
}

AbaqusCDPStateIntegrator::Result
AbaqusCDPStateIntegrator::assembleBackboneResult(
    const SymmetricTensor & total_strain,
    const double time_step,
    const State & old_state,
    const AbaqusCDPLocalIntegrator::Result & backbone) const
{
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::STATE);
  validateAssemblyInputs(total_strain, time_step, old_state);
  return assembleResult(total_strain, time_step, old_state, backbone);
}

AbaqusCDPStateIntegrator::SymmetricTensor
AbaqusCDPStateIntegrator::backboneEffectiveStress(const SymmetricTensor & total_strain,
                                                  const State & state) const
{
  return _backbone_integrator.elasticStress(
      subtractStateTensor(total_strain, state.backbone.plastic_strain));
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

  const auto viscous_effective_stress = _backbone_integrator.elasticStress(
      subtractStateTensor(total_strain, new_state.viscous_plastic_strain));
  const double tension_weight =
      AbaqusCDPFormula::stressInvariants(viscous_effective_stress).tension_weight;
  const auto damage = AbaqusCDPFormula::combineDamage(new_state.viscous_compression_damage,
                                                      new_state.viscous_tension_damage,
                                                      _parameters.tension_recovery,
                                                      _parameters.compression_recovery,
                                                      tension_weight);
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
          backbone.backbone_compression_damage - new_state.viscous_compression_damage};
}

AbaqusCDPStateIntegrator::LinearizedResult
AbaqusCDPStateIntegrator::integrateLinearized(const SymmetricTensor & total_strain,
                                              const double time_step,
                                              const State & old_state) const
{
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::STATE_LINEARIZED);
  validateAssemblyInputs(total_strain, time_step, old_state);
  const auto backbone = integrateBackboneLinearized(total_strain, old_state.backbone);
  return assembleBackboneLinearized(total_strain, time_step, old_state, backbone);
}

AbaqusCDPStateIntegrator::LinearizedResult
AbaqusCDPStateIntegrator::assembleBackboneLinearized(
    const SymmetricTensor & total_strain,
    const double time_step,
    const State & old_state,
    const AbaqusCDPLocalIntegrator::LinearizedResult & backbone) const
{
  validateAssemblyInputs(total_strain, time_step, old_state);
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

  SymmetricTensor tension_weight_gradient;
  double stress_scale = 1.0;
  for (const double value : linearized.result.viscous_effective_stress)
    stress_scale = std::max(stress_scale, std::abs(value));
  const double weight_step = 1.0e-7 * stress_scale;
  for (std::size_t column = 0; column < 6; ++column)
  {
    auto plus = linearized.result.viscous_effective_stress;
    auto minus = linearized.result.viscous_effective_stress;
    plus[column] += weight_step;
    minus[column] -= weight_step;
    tension_weight_gradient[column] =
        (AbaqusCDPFormula::stressInvariants(plus).tension_weight -
         AbaqusCDPFormula::stressInvariants(minus).tension_weight) /
        (2.0 * weight_step);
  }

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
    if (input >= 6)
    {
      const std::size_t old_state_index = input - 6;
      if (old_state_index >= 8 && old_state_index < 14)
        old_viscous_plastic_derivative[old_state_index - 8] = 1.0;
      else if (old_state_index == 14)
        old_tension_damage_derivative = 1.0;
      else if (old_state_index == 15)
        old_compression_damage_derivative = 1.0;
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
    for (std::size_t row = 0; row < 6; ++row)
      tension_weight_derivative +=
          tension_weight_gradient[row] * effective_stress_derivative[row];

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
    const double stiffness_derivative =
        compression_factor_derivative * tension_factor +
        compression_factor * tension_factor_derivative;

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
  }
  return linearized;
}
