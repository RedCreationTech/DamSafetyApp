#include "CDPDiagnostics.h"
#include "AbaqusCDPSubstepIntegrator.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
using SymmetricTensor = AbaqusCDPSubstepIntegrator::SymmetricTensor;

[[noreturn]] void
substepError(const std::string & message)
{
  throw std::runtime_error("AbaqusCDPSubstepIntegrator: " + message);
}

bool
finiteSubstepTensor(const SymmetricTensor & tensor)
{
  return std::all_of(
      tensor.begin(), tensor.end(), [](const double value) { return std::isfinite(value); });
}

SymmetricTensor
substepDifference(const SymmetricTensor & left, const SymmetricTensor & right)
{
  SymmetricTensor result;
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = left[i] - right[i];
  return result;
}

SymmetricTensor
substepInterpolate(const SymmetricTensor & start,
                   const SymmetricTensor & increment,
                   const double fraction)
{
  SymmetricTensor result;
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = start[i] + fraction * increment[i];
  return result;
}

double
maximumAbsoluteComponent(const SymmetricTensor & tensor)
{
  double result = 0.0;
  for (const double value : tensor)
    result = std::max(result, std::abs(value));
  return result;
}

bool
isPowerOfTwo(const unsigned int value)
{
  return value != 0 && (value & (value - 1)) == 0;
}
}

AbaqusCDPSubstepIntegrator::AbaqusCDPSubstepIntegrator(
    const AbaqusCDPStateIntegrator & state_integrator,
    Parameters parameters)
  : _state_integrator(state_integrator), _parameters(parameters)
{
  if (!isPowerOfTwo(_parameters.maximum_substeps))
    substepError("maximum_substeps must be a positive power of two");
  if (!std::isfinite(_parameters.maximum_strain_increment) ||
      _parameters.maximum_strain_increment < 0.0)
    substepError("maximum_strain_increment must be finite and nonnegative");
  if (!std::isfinite(_parameters.tangent_perturbation) ||
      _parameters.tangent_perturbation <= 0.0)
    substepError("tangent_perturbation must be finite and positive");
}

AbaqusCDPSubstepIntegrator::Result
AbaqusCDPSubstepIntegrator::integrate(const SymmetricTensor & old_total_strain,
                                      const SymmetricTensor & new_total_strain,
                                      const double time_step,
                                      const State & old_state) const
{
  if (!finiteSubstepTensor(old_total_strain) || !finiteSubstepTensor(new_total_strain))
    substepError("old or new total strain contains a non-finite value");
  if (!std::isfinite(time_step) || time_step < 0.0)
    substepError("time step must be finite and nonnegative");

  const auto total_increment = substepDifference(new_total_strain, old_total_strain);
  unsigned int substeps = 1;
  if (_parameters.maximum_strain_increment > 0.0)
    while (maximumAbsoluteComponent(total_increment) / substeps >
           _parameters.maximum_strain_increment)
    {
      if (substeps == _parameters.maximum_substeps)
        substepError("proactive strain-increment limit requires more than maximum_substeps");
      substeps *= 2;
    }
  const bool proactively_partitioned = substeps > 1;

  unsigned int cutback_count = 0;
  unsigned int attempted_partitions = 0;
  std::string last_error;
  unsigned int last_failed_substep = 0;
  unsigned int last_partition = substeps;

  while (substeps <= _parameters.maximum_substeps)
  {
    CDPDiagnostics::Scope partition_scope(CDPDiagnostics::PARTITION);
    if (CDPDiagnostics::current) CDPDiagnostics::current->partition=substeps;
    ++attempted_partitions;
    State working_state = old_state;
    std::optional<AbaqusCDPStateIntegrator::Result> final_result;
    unsigned int total_local_iterations = 0;
    unsigned int total_jacobian_fallbacks = 0;
    unsigned int total_automatic_jacobian_evaluations = 0;
    unsigned int total_finite_difference_jacobian_evaluations = 0;
    unsigned int total_local_factorizations = 0;
    unsigned int total_local_backsolves = 0;
    bool partition_succeeded = true;

    for (unsigned int i = 1; i <= substeps; ++i)
    {
      const auto target = substepInterpolate(
          old_total_strain, total_increment, static_cast<double>(i) / substeps);
      if (CDPDiagnostics::current) CDPDiagnostics::current->substep=i;
      try
      {
        auto step_result = _state_integrator.integrate(
            target, time_step / static_cast<double>(substeps), working_state);
        total_local_iterations += step_result.backbone.iterations;
        total_jacobian_fallbacks += step_result.backbone.jacobian_fallbacks;
        total_automatic_jacobian_evaluations +=
            step_result.backbone.automatic_jacobian_evaluations;
        total_finite_difference_jacobian_evaluations +=
            step_result.backbone.finite_difference_jacobian_evaluations;
        total_local_factorizations += step_result.backbone.local_factorizations;
        total_local_backsolves += step_result.backbone.local_backsolves;
        working_state = step_result.state;
        final_result = std::move(step_result);
      }
      catch (const std::exception & error)
      {
        partition_succeeded = false;
        last_error = error.what();
        last_failed_substep = i;
        last_partition = substeps;
        break;
      }
    }

    if (partition_succeeded && final_result)
      return {std::move(*final_result),
              substeps,
              cutback_count,
              attempted_partitions,
              total_local_iterations,
              total_jacobian_fallbacks,
              total_automatic_jacobian_evaluations,
              total_finite_difference_jacobian_evaluations,
              total_local_factorizations,
              total_local_backsolves,
              proactively_partitioned};

    partition_scope.failed();
    if (substeps == _parameters.maximum_substeps)
      break;
    substeps *= 2;
    ++cutback_count;
  }

  std::ostringstream message;
  message << "failed after partition=" << last_partition
          << ", failed_substep=" << last_failed_substep
          << ", attempted_partitions=" << attempted_partitions
          << ", last_error=" << last_error;
  substepError(message.str());
}

AbaqusCDPSubstepIntegrator::LinearizedResult
AbaqusCDPSubstepIntegrator::integrateLinearized(const SymmetricTensor & old_total_strain,
                                                const SymmetricTensor & new_total_strain,
                                                const double time_step,
                                                const State & old_state) const
{
  if (!finiteSubstepTensor(old_total_strain) || !finiteSubstepTensor(new_total_strain))
    substepError("old or new total strain contains a non-finite value");
  if (!std::isfinite(time_step) || time_step < 0.0)
    substepError("time step must be finite and nonnegative");

  const auto total_increment = substepDifference(new_total_strain, old_total_strain);
  unsigned int substeps = 1;
  if (_parameters.maximum_strain_increment > 0.0)
    while (maximumAbsoluteComponent(total_increment) / substeps >
           _parameters.maximum_strain_increment)
    {
      if (substeps == _parameters.maximum_substeps)
        substepError("proactive strain-increment limit requires more than maximum_substeps");
      substeps *= 2;
    }
  const bool proactively_partitioned = substeps > 1;

  unsigned int cutback_count = 0;
  unsigned int attempted_partitions = 0;
  std::string last_error;
  unsigned int last_failed_substep = 0;
  unsigned int last_partition = substeps;

  while (substeps <= _parameters.maximum_substeps)
  {
    CDPDiagnostics::Scope partition_scope(CDPDiagnostics::PARTITION);
    if (CDPDiagnostics::current) CDPDiagnostics::current->partition=substeps;
    ++attempted_partitions;
    State working_state = old_state;
    std::optional<AbaqusCDPStateIntegrator::LinearizedResult> final_result;
    std::array<std::array<double, AbaqusCDPStateIntegrator::state_size>, 6>
        state_sensitivity = {};
    TangentMatrix tangent = {};
    unsigned int total_local_iterations = 0;
    unsigned int total_jacobian_fallbacks = 0;
    unsigned int total_automatic_jacobian_evaluations = 0;
    unsigned int total_finite_difference_jacobian_evaluations = 0;
    unsigned int total_local_factorizations = 0;
    unsigned int total_local_backsolves = 0;
    bool partition_succeeded = true;

    for (unsigned int i = 1; i <= substeps; ++i)
    {
      const double fraction = static_cast<double>(i) / substeps;
      const auto target = substepInterpolate(old_total_strain, total_increment, fraction);
      if (CDPDiagnostics::current) CDPDiagnostics::current->substep=i;
      try
      {
        auto step_result = _state_integrator.integrateLinearized(
            target, time_step / static_cast<double>(substeps), working_state);
        total_local_iterations += step_result.result.backbone.iterations;
        total_jacobian_fallbacks += step_result.result.backbone.jacobian_fallbacks;
        total_automatic_jacobian_evaluations +=
            step_result.result.backbone.automatic_jacobian_evaluations;
        total_finite_difference_jacobian_evaluations +=
            step_result.result.backbone.finite_difference_jacobian_evaluations;
        total_local_factorizations += step_result.result.backbone.local_factorizations;
        total_local_backsolves += step_result.result.backbone.local_backsolves;

        std::array<std::array<double, AbaqusCDPStateIntegrator::state_size>, 6>
            new_state_sensitivity = {};
        TangentMatrix new_tangent = {};
        {
        CDPDiagnostics::Scope chain_scope(CDPDiagnostics::STATE_CHAIN);
        for (std::size_t final_column = 0; final_column < 6; ++final_column)
        {
          std::array<double, AbaqusCDPStateIntegrator::transition_size> input_derivative = {};
          input_derivative[final_column] = fraction;
          for (std::size_t state_row = 0;
               state_row < AbaqusCDPStateIntegrator::state_size;
               ++state_row)
            input_derivative[6 + state_row] = state_sensitivity[final_column][state_row];

          for (std::size_t output_row = 0; output_row < 6; ++output_row)
            for (std::size_t input = 0;
                 input < AbaqusCDPStateIntegrator::transition_size;
                 ++input)
              new_tangent[final_column][output_row] +=
                  step_result.derivative[input][output_row] * input_derivative[input];
          for (std::size_t state_row = 0;
               state_row < AbaqusCDPStateIntegrator::state_size;
               ++state_row)
            for (std::size_t input = 0;
                 input < AbaqusCDPStateIntegrator::transition_size;
                 ++input)
              new_state_sensitivity[final_column][state_row] +=
                  step_result.derivative[input][6 + state_row] * input_derivative[input];
        }
        }
        if (CDPDiagnostics::current && CDPDiagnostics::current->trace)
        {
          const double start_fraction = static_cast<double>(i - 1) / substeps;
          const auto start_target =
              substepInterpolate(old_total_strain, total_increment, start_fraction);
          const auto start_backbone_stress =
              _state_integrator.backboneEffectiveStress(start_target, working_state);
          const auto plastic_increment = substepDifference(
              step_result.result.state.backbone.plastic_strain,
              working_state.backbone.plastic_strain);
          const double delta_kappa_t =
              step_result.result.state.backbone.tensile_equivalent_plastic_strain -
              working_state.backbone.tensile_equivalent_plastic_strain;
          const double delta_kappa_c =
              step_result.result.state.backbone.compressive_equivalent_plastic_strain -
              working_state.backbone.compressive_equivalent_plastic_strain;
          const auto start_invariants = AbaqusCDPFormula::stressInvariants(start_backbone_stress);
          const auto end_invariants =
              AbaqusCDPFormula::stressInvariants(step_result.result.backbone.effective_stress);
          const auto plastic_increment_principal =
              AbaqusCDPFormula::stressInvariants(plastic_increment).principal_stress;
          const double positive_plastic_principal_increment =
              std::max(0.0, plastic_increment_principal[2]);
          const double start_backbone_tension_weight = start_invariants.tension_weight;
          const double backbone_tension_weight = end_invariants.tension_weight;
          const bool maximum_principal_zero_crossing =
              (start_invariants.maximum_principal_stress < 0.0 &&
               end_invariants.maximum_principal_stress >= 0.0) ||
              (start_invariants.maximum_principal_stress > 0.0 &&
               end_invariants.maximum_principal_stress <= 0.0);
          double maximum_principal_zero_crossing_fraction = -1.0;
          if (maximum_principal_zero_crossing)
          {
            const double denominator = start_invariants.maximum_principal_stress -
                                       end_invariants.maximum_principal_stress;
            if (std::abs(denominator) > 0.0)
              maximum_principal_zero_crossing_fraction =
                  start_invariants.maximum_principal_stress / denominator;
          }
          const double viscous_tension_weight =
              AbaqusCDPFormula::stressInvariants(step_result.result.viscous_effective_stress)
                  .tension_weight;
          std::ostringstream payload;
          payload<<"\"target\":";CDPDiagnostics::json(payload,target);
          payload<<",\"start_target\":";CDPDiagnostics::json(payload,start_target);
          payload<<",\"old_state\":";CDPDiagnostics::stateJson(payload,working_state);
          payload<<",\"new_state\":";CDPDiagnostics::stateJson(payload,step_result.result.state);
          payload<<",\"history_integration_rule\":\"end\"";
          payload<<",\"start_backbone_stress\":";CDPDiagnostics::json(payload,start_backbone_stress);
          payload<<",\"backbone_stress\":";CDPDiagnostics::json(payload,step_result.result.backbone.effective_stress);
          payload<<",\"start_principal_stress\":";CDPDiagnostics::json(payload,start_invariants.principal_stress);
          payload<<",\"end_principal_stress\":";CDPDiagnostics::json(payload,end_invariants.principal_stress);
          payload<<",\"viscous_stress\":";CDPDiagnostics::json(payload,step_result.result.viscous_effective_stress);
          payload<<",\"plastic_increment\":";CDPDiagnostics::json(payload,plastic_increment);
          payload<<",\"plastic_increment_principal\":";CDPDiagnostics::json(payload,plastic_increment_principal);
          payload<<",\"positive_plastic_principal_increment\":";CDPDiagnostics::json(payload,positive_plastic_principal_increment);
          payload<<",\"delta_kappa_t\":";CDPDiagnostics::json(payload,delta_kappa_t);
          payload<<",\"delta_kappa_t_end_rule\":";CDPDiagnostics::json(payload,delta_kappa_t);
          payload<<",\"delta_kappa_c\":";CDPDiagnostics::json(payload,delta_kappa_c);
          payload<<",\"start_backbone_tension_weight\":";CDPDiagnostics::json(payload,start_backbone_tension_weight);
          payload<<",\"backbone_tension_weight\":";CDPDiagnostics::json(payload,backbone_tension_weight);
          payload<<",\"end_backbone_tension_weight\":";CDPDiagnostics::json(payload,backbone_tension_weight);
          payload<<",\"maximum_principal_zero_crossing\":"<<(maximum_principal_zero_crossing ? "true" : "false");
          payload<<",\"maximum_principal_zero_crossing_fraction\":";CDPDiagnostics::json(payload,maximum_principal_zero_crossing_fraction);
          payload<<",\"viscous_tension_weight\":";CDPDiagnostics::json(payload,viscous_tension_weight);
          payload<<",\"active_branch\":\""<<AbaqusCDPLocalIntegrator::branchName(step_result.result.backbone.active_branch)<<"\"";
          payload<<",\"backbone_damage_t\":";CDPDiagnostics::json(payload,step_result.result.backbone.backbone_tension_damage);
          payload<<",\"backbone_damage_c\":";CDPDiagnostics::json(payload,step_result.result.backbone.backbone_compression_damage);
          payload<<",\"viscous_damage_t\":";CDPDiagnostics::json(payload,step_result.result.state.viscous_tension_damage);
          payload<<",\"viscous_damage_c\":";CDPDiagnostics::json(payload,step_result.result.state.viscous_compression_damage);
          payload<<",\"combined_damage\":";CDPDiagnostics::json(payload,step_result.result.damage.damage);
          payload<<",\"plastic_strain_lag_norm\":";CDPDiagnostics::json(payload,step_result.result.plastic_strain_lag_norm);
          payload<<",\"tension_damage_lag\":";CDPDiagnostics::json(payload,step_result.result.tension_damage_lag);
          payload<<",\"compression_damage_lag\":";CDPDiagnostics::json(payload,step_result.result.compression_damage_lag);
          payload<<",\"dt_over_relaxation_time\":";CDPDiagnostics::json(payload,step_result.result.dt_over_relaxation_time);
          payload<<",\"yield\":";CDPDiagnostics::json(payload,step_result.result.backbone.final_yield);
          payload<<",\"residual\":";CDPDiagnostics::json(payload,step_result.result.backbone.residual_norm);
          payload<<",\"multiplier\":";CDPDiagnostics::json(payload,step_result.result.backbone.plastic_multiplier);
          CDPDiagnostics::event("substep",payload.str());
        }
        state_sensitivity = new_state_sensitivity;
        tangent = new_tangent;
        working_state = step_result.result.state;
        final_result = std::move(step_result);
      }
      catch (const std::exception & error)
      {
        partition_succeeded = false;
        last_error = error.what();
        last_failed_substep = i;
        last_partition = substeps;
        break;
      }
    }

    if (partition_succeeded && final_result)
    {
      Result result{std::move(final_result->result),
                    substeps,
                    cutback_count,
                    attempted_partitions,
                    total_local_iterations,
                    total_jacobian_fallbacks,
                    total_automatic_jacobian_evaluations,
                    total_finite_difference_jacobian_evaluations,
                    total_local_factorizations,
                    total_local_backsolves,
                    proactively_partitioned};
      return {std::move(result), tangent};
    }

    partition_scope.failed();
    if (substeps == _parameters.maximum_substeps)
      break;
    substeps *= 2;
    ++cutback_count;
  }

  std::ostringstream message;
  message << "linearization failed after partition=" << last_partition
          << ", failed_substep=" << last_failed_substep
          << ", attempted_partitions=" << attempted_partitions
          << ", last_error=" << last_error;
  substepError(message.str());
}

AbaqusCDPSubstepIntegrator::ReferenceTangent
AbaqusCDPSubstepIntegrator::referenceTangent(const SymmetricTensor & old_total_strain,
                                             const SymmetricTensor & new_total_strain,
                                             const double time_step,
                                             const State & old_state) const
{
  ReferenceTangent result;
  result.perturbation = _parameters.tangent_perturbation;
  for (std::size_t column = 0; column < 6; ++column)
  {
    auto plus_strain = new_total_strain;
    auto minus_strain = new_total_strain;
    plus_strain[column] += result.perturbation;
    minus_strain[column] -= result.perturbation;
    const auto plus = integrate(old_total_strain, plus_strain, time_step, old_state);
    const auto minus = integrate(old_total_strain, minus_strain, time_step, old_state);
    result.plus_substeps[column] = plus.accepted_substeps;
    result.minus_substeps[column] = minus.accepted_substeps;
    for (std::size_t row = 0; row < 6; ++row)
      result.value[column][row] =
          (plus.final_result.cauchy_stress[row] - minus.final_result.cauchy_stress[row]) /
          (2.0 * result.perturbation);
  }
  return result;
}

AbaqusCDPSubstepIntegrator::SymmetricTensor
AbaqusCDPSubstepIntegrator::directionalDerivative(const SymmetricTensor & old_total_strain,
                                                   const SymmetricTensor & new_total_strain,
                                                   const double time_step,
                                                   const State & old_state,
                                                   const SymmetricTensor & direction,
                                                   const double perturbation) const
{
  if (!finiteSubstepTensor(direction) || !std::isfinite(perturbation) || perturbation <= 0.0)
    substepError("direction and directional perturbation must be finite and positive");
  auto plus_strain = new_total_strain;
  auto minus_strain = new_total_strain;
  for (std::size_t i = 0; i < 6; ++i)
  {
    plus_strain[i] += perturbation * direction[i];
    minus_strain[i] -= perturbation * direction[i];
  }
  const auto plus = integrate(old_total_strain, plus_strain, time_step, old_state);
  const auto minus = integrate(old_total_strain, minus_strain, time_step, old_state);
  SymmetricTensor result;
  for (std::size_t i = 0; i < 6; ++i)
    result[i] = (plus.final_result.cauchy_stress[i] - minus.final_result.cauchy_stress[i]) /
                (2.0 * perturbation);
  return result;
}

AbaqusCDPSubstepIntegrator::SymmetricTensor
AbaqusCDPSubstepIntegrator::applyTangent(const TangentMatrix & tangent,
                                         const SymmetricTensor & direction)
{
  if (!finiteSubstepTensor(direction))
    substepError("tangent direction contains a non-finite value");
  SymmetricTensor result = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (std::size_t column = 0; column < 6; ++column)
    for (std::size_t row = 0; row < 6; ++row)
      result[row] += tangent[column][row] * direction[column];
  return result;
}
