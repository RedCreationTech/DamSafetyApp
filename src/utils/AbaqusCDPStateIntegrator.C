#include "AbaqusCDPStateIntegrator.h"

#include <algorithm>
#include <cmath>
#include <limits>
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
  if (!finiteStateTensor(total_strain) || !finiteStateTensor(old_state.viscous_plastic_strain))
    stateIntegrationError("total strain or old viscous plastic strain contains a non-finite value");
  if (!std::isfinite(time_step) || time_step < 0.0)
    stateIntegrationError("time step must be finite and nonnegative");
  validateDamage(old_state.viscous_tension_damage, "old viscous tension damage");
  validateDamage(old_state.viscous_compression_damage, "old viscous compression damage");

  // The old state remains immutable. Any failure below leaves the caller's
  // checkpoint untouched and makes a cutback retry deterministic.
  const auto backbone = _backbone_integrator.integrate(total_strain, old_state.backbone);
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
  const double dt_over_relaxation_time =
      _parameters.relaxation_time == 0.0
          ? std::numeric_limits<double>::infinity()
          : time_step / _parameters.relaxation_time;

  return {backbone,
          viscous_effective_stress,
          cauchy_stress,
          new_state,
          damage,
          dt_over_relaxation_time,
          stateTensorNorm(plastic_lag),
          backbone.backbone_tension_damage - new_state.viscous_tension_damage,
          backbone.backbone_compression_damage - new_state.viscous_compression_damage};
}
