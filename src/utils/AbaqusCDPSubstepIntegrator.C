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

double
positivePrincipalMeasure(const SymmetricTensor & tensor)
{
  return std::max(0.0, AbaqusCDPFormula::stressInvariants(tensor).principal_stress.back());
}

double
negativePrincipalMeasure(const SymmetricTensor & tensor)
{
  return std::max(0.0, -AbaqusCDPFormula::stressInvariants(tensor).principal_stress.front());
}

template <typename Function>
SymmetricTensor
scalarTensorGradient(const SymmetricTensor & argument, const double step, Function function)
{
  SymmetricTensor gradient = {};
  for (std::size_t component = 0; component < gradient.size(); ++component)
  {
    auto plus = argument;
    auto minus = argument;
    plus[component] += step;
    minus[component] -= step;
    gradient[component] = (function(plus) - function(minus)) / (2.0 * step);
  }
  return gradient;
}

double
tensorDot(const SymmetricTensor & left, const SymmetricTensor & right)
{
  double value = 0.0;
  for (std::size_t component = 0; component < left.size(); ++component)
    value += left[component] * right[component];
  return value;
}

bool
isPowerOfTwo(const unsigned int value)
{
  return value != 0 && (value & (value - 1)) == 0;
}

struct HistoryPathAccumulator
{
  double tensile_measure = 0.0;
  double tensile_weighted_measure = 0.0;
  double compressive_measure = 0.0;
  double compressive_weighted_measure = 0.0;

  void add(const SymmetricTensor & old_plastic_strain,
           const SymmetricTensor & new_plastic_strain,
           const SymmetricTensor & start_effective_stress,
           const SymmetricTensor & end_effective_stress)
  {
    const auto plastic_increment = substepDifference(new_plastic_strain, old_plastic_strain);
    const auto plastic_principal =
        AbaqusCDPFormula::stressInvariants(plastic_increment).principal_stress;
    const double positive_measure = std::max(0.0, plastic_principal.back());
    const double negative_measure = std::max(0.0, -plastic_principal.front());
    const double start_weight =
        AbaqusCDPFormula::stressInvariants(start_effective_stress).tension_weight;
    const double end_weight =
        AbaqusCDPFormula::stressInvariants(end_effective_stress).tension_weight;
    const double average_weight = 0.5 * (start_weight + end_weight);

    tensile_measure += positive_measure;
    tensile_weighted_measure += average_weight * positive_measure;
    compressive_measure += negative_measure;
    compressive_weighted_measure += (1.0 - average_weight) * negative_measure;
  }

  double tensileWeight(const double fallback) const
  {
    return tensile_measure > 0.0 ? tensile_weighted_measure / tensile_measure : fallback;
  }

  double compressiveWeight(const double fallback) const
  {
    return compressive_measure > 0.0 ? compressive_weighted_measure / compressive_measure
                                     : fallback;
  }
};
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
  if (_parameters.aggregate_backbone_history_to_global_step &&
      !_parameters.defer_viscous_update_to_global_step)
    substepError("aggregate_backbone_history_to_global_step requires "
                 "defer_viscous_update_to_global_step");
  if (_parameters.integrate_backbone_history_weights_over_substeps &&
      !_parameters.aggregate_backbone_history_to_global_step)
    substepError("integrate_backbone_history_weights_over_substeps requires "
                 "aggregate_backbone_history_to_global_step");
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
    std::optional<AbaqusCDPLocalIntegrator::Result> final_backbone;
    HistoryPathAccumulator history_path;
    SymmetricTensor previous_effective_stress =
        _state_integrator.backboneEffectiveStress(old_total_strain, old_state);
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
        if (_parameters.defer_viscous_update_to_global_step)
        {
          const auto old_backbone = working_state.backbone;
          auto backbone =
              _state_integrator.integrateBackbone(target, working_state.backbone);
          total_local_iterations += backbone.iterations;
          total_jacobian_fallbacks += backbone.jacobian_fallbacks;
          total_automatic_jacobian_evaluations +=
              backbone.automatic_jacobian_evaluations;
          total_finite_difference_jacobian_evaluations +=
              backbone.finite_difference_jacobian_evaluations;
          total_local_factorizations += backbone.local_factorizations;
          total_local_backsolves += backbone.local_backsolves;
          if (_parameters.integrate_backbone_history_weights_over_substeps)
            history_path.add(old_backbone.plastic_strain,
                             backbone.state.plastic_strain,
                             previous_effective_stress,
                             backbone.effective_stress);
          previous_effective_stress = backbone.effective_stress;
          working_state.backbone = backbone.state;
          final_backbone = std::move(backbone);
        }
        else
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

    if (partition_succeeded && _parameters.defer_viscous_update_to_global_step && final_backbone)
      try
      {
        if (_parameters.aggregate_backbone_history_to_global_step)
        {
          if (_parameters.integrate_backbone_history_weights_over_substeps)
          {
            const double final_tension_weight =
                AbaqusCDPFormula::stressInvariants(final_backbone->effective_stress).tension_weight;
            *final_backbone = _state_integrator.aggregateBackboneHistory(
                old_state.backbone,
                *final_backbone,
                history_path.tensileWeight(final_tension_weight),
                history_path.compressiveWeight(1.0 - final_tension_weight));
          }
          else
            *final_backbone = _state_integrator.aggregateBackboneHistory(
                old_state.backbone, *final_backbone);
        }
        final_result = _state_integrator.assembleBackboneResult(
            new_total_strain, time_step, old_state, *final_backbone);
      }
      catch (const std::exception & error)
      {
        partition_succeeded = false;
        last_error = error.what();
        last_failed_substep = substeps;
        last_partition = substeps;
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
AbaqusCDPSubstepIntegrator::integrateDeferredViscousLinearized(
    const SymmetricTensor & old_total_strain,
    const SymmetricTensor & new_total_strain,
    const double time_step,
    const State & old_state) const
{
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::STATE_LINEARIZED);
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
    if (CDPDiagnostics::current)
      CDPDiagnostics::current->partition = substeps;
    ++attempted_partitions;

    auto working_backbone = old_state.backbone;
    std::optional<AbaqusCDPLocalIntegrator::LinearizedResult> final_backbone;
    AbaqusCDPLocalIntegrator::TransitionJacobian chained_derivative = {};
    std::array<std::array<double, 8>, AbaqusCDPLocalIntegrator::transition_size>
        backbone_state_sensitivity = {};
    for (std::size_t state = 0; state < 8; ++state)
      backbone_state_sensitivity[6 + state][state] = 1.0;

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
      if (CDPDiagnostics::current)
        CDPDiagnostics::current->substep = i;
      try
      {
        const auto old_backbone = working_backbone;
        auto step = _state_integrator.integrateBackboneLinearized(target, working_backbone);
        total_local_iterations += step.result.iterations;
        total_jacobian_fallbacks += step.result.jacobian_fallbacks;
        total_automatic_jacobian_evaluations +=
            step.result.automatic_jacobian_evaluations;
        total_finite_difference_jacobian_evaluations +=
            step.result.finite_difference_jacobian_evaluations;
        total_local_factorizations += step.result.local_factorizations;
        total_local_backsolves += step.result.local_backsolves;

        std::array<std::array<double, 8>, AbaqusCDPLocalIntegrator::transition_size>
            new_state_sensitivity = {};
        AbaqusCDPLocalIntegrator::TransitionJacobian new_chained_derivative = {};
        {
          CDPDiagnostics::Scope chain_scope(CDPDiagnostics::STATE_CHAIN);
          for (std::size_t original_input = 0;
               original_input < AbaqusCDPLocalIntegrator::transition_size;
               ++original_input)
          {
            std::array<double, AbaqusCDPLocalIntegrator::transition_size> input_derivative = {};
            if (original_input < 6)
              input_derivative[original_input] = fraction;
            for (std::size_t state = 0; state < 8; ++state)
              input_derivative[6 + state] =
                  backbone_state_sensitivity[original_input][state];

            for (std::size_t output = 0;
                 output < AbaqusCDPLocalIntegrator::transition_size;
                 ++output)
              for (std::size_t input = 0;
                   input < AbaqusCDPLocalIntegrator::transition_size;
                   ++input)
                new_chained_derivative[original_input][output] +=
                    step.derivative[input][output] * input_derivative[input];
            for (std::size_t state = 0; state < 8; ++state)
              new_state_sensitivity[original_input][state] =
                  new_chained_derivative[original_input][6 + state];
          }
        }

        if (CDPDiagnostics::current && CDPDiagnostics::current->trace)
        {
          const double start_fraction = static_cast<double>(i - 1) / substeps;
          const auto start_target =
              substepInterpolate(old_total_strain, total_increment, start_fraction);
          auto start_state = old_state;
          start_state.backbone = old_backbone;
          const auto start_effective_stress =
              _state_integrator.backboneEffectiveStress(start_target, start_state);
          const auto plastic_increment =
              substepDifference(step.result.state.plastic_strain, old_backbone.plastic_strain);
          const auto plastic_increment_invariants =
              AbaqusCDPFormula::stressInvariants(plastic_increment);
          const auto start_invariants =
              AbaqusCDPFormula::stressInvariants(start_effective_stress);
          const auto end_invariants =
              AbaqusCDPFormula::stressInvariants(step.result.effective_stress);
          const double positive_plastic_measure =
              std::max(0.0, plastic_increment_invariants.principal_stress.back());
          const double negative_plastic_measure =
              std::max(0.0, -plastic_increment_invariants.principal_stress.front());
          const double raw_delta_kappa_t =
              end_invariants.tension_weight * positive_plastic_measure;
          const double raw_delta_kappa_c =
              (1.0 - end_invariants.tension_weight) * negative_plastic_measure;
          SymmetricTensor flow_direction = {};
          if (step.result.plastic_multiplier > 0.0)
            for (std::size_t component = 0; component < flow_direction.size(); ++component)
              flow_direction[component] =
                  plastic_increment[component] / step.result.plastic_multiplier;
          const auto flow_direction_principal =
              AbaqusCDPFormula::stressInvariants(flow_direction).principal_stress;
          std::ostringstream payload;
          payload << "\"target\":";
          CDPDiagnostics::json(payload, target);
          payload << ",\"start_target\":";
          CDPDiagnostics::json(payload, start_target);
          payload << ",\"history_integration_rule\":\"material_substep_end\"";
          payload << ",\"viscous_state_commit\":\"deferred_to_global_accepted_step\"";
          payload << ",\"start_effective_stress\":";
          CDPDiagnostics::json(payload, start_effective_stress);
          payload << ",\"end_effective_stress\":";
          CDPDiagnostics::json(payload, step.result.effective_stress);
          payload << ",\"start_effective_principal\":";
          CDPDiagnostics::json(payload, start_invariants.principal_stress);
          payload << ",\"end_effective_principal\":";
          CDPDiagnostics::json(payload, end_invariants.principal_stress);
          payload << ",\"start_tension_weight\":";
          CDPDiagnostics::json(payload, start_invariants.tension_weight);
          payload << ",\"end_tension_weight\":";
          CDPDiagnostics::json(payload, end_invariants.tension_weight);
          payload << ",\"plastic_increment\":";
          CDPDiagnostics::json(payload, plastic_increment);
          payload << ",\"plastic_increment_principal\":";
          CDPDiagnostics::json(payload, plastic_increment_invariants.principal_stress);
          payload << ",\"positive_plastic_measure\":";
          CDPDiagnostics::json(payload, positive_plastic_measure);
          payload << ",\"negative_plastic_measure\":";
          CDPDiagnostics::json(payload, negative_plastic_measure);
          payload << ",\"flow_direction\":";
          CDPDiagnostics::json(payload, flow_direction);
          payload << ",\"flow_direction_principal\":";
          CDPDiagnostics::json(payload, flow_direction_principal);
          payload << ",\"raw_delta_kappa_t\":";
          CDPDiagnostics::json(payload, raw_delta_kappa_t);
          payload << ",\"raw_delta_kappa_c\":";
          CDPDiagnostics::json(payload, raw_delta_kappa_c);
          payload << ",\"old_kappa_t\":";
          CDPDiagnostics::json(payload, old_backbone.tensile_equivalent_plastic_strain);
          payload << ",\"new_kappa_t\":";
          CDPDiagnostics::json(payload,
                               step.result.state.tensile_equivalent_plastic_strain);
          payload << ",\"old_kappa_c\":";
          CDPDiagnostics::json(payload, old_backbone.compressive_equivalent_plastic_strain);
          payload << ",\"new_kappa_c\":";
          CDPDiagnostics::json(payload,
                               step.result.state.compressive_equivalent_plastic_strain);
          payload << ",\"backbone_damage_t\":";
          CDPDiagnostics::json(payload, step.result.backbone_tension_damage);
          payload << ",\"backbone_damage_c\":";
          CDPDiagnostics::json(payload, step.result.backbone_compression_damage);
          payload << ",\"active_branch\":\""
                  << AbaqusCDPLocalIntegrator::branchName(step.result.active_branch) << "\"";
          payload << ",\"trial_yield\":";
          CDPDiagnostics::json(payload, step.result.trial_yield);
          payload << ",\"final_yield\":";
          CDPDiagnostics::json(payload, step.result.final_yield);
          payload << ",\"residual_norm\":";
          CDPDiagnostics::json(payload, step.result.residual_norm);
          payload << ",\"plastic_multiplier\":";
          CDPDiagnostics::json(payload, step.result.plastic_multiplier);
          CDPDiagnostics::event("substep", payload.str());
        }

        working_backbone = step.result.state;
        backbone_state_sensitivity = new_state_sensitivity;
        chained_derivative = new_chained_derivative;
        final_backbone = std::move(step);
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

    if (partition_succeeded && final_backbone)
    {
      try
      {
        final_backbone->derivative = chained_derivative;
        if (_parameters.aggregate_backbone_history_to_global_step)
        {
          const auto plastic_increment = substepDifference(
              final_backbone->result.state.plastic_strain, old_state.backbone.plastic_strain);
          const double tension_weight =
              AbaqusCDPFormula::stressInvariants(final_backbone->result.effective_stress)
                  .tension_weight;
          const double tensile_measure = positivePrincipalMeasure(plastic_increment);
          const double compressive_measure = negativePrincipalMeasure(plastic_increment);

          double stress_scale = 1.0;
          for (const double value : final_backbone->result.effective_stress)
            stress_scale = std::max(stress_scale, std::abs(value));
          const double stress_step = 1.0e-7 * stress_scale;
          const double plastic_step = std::max(
              _parameters.tangent_perturbation,
              1.0e-7 * std::max(maximumAbsoluteComponent(plastic_increment), 1.0e-12));
          const auto tension_weight_gradient = scalarTensorGradient(
              final_backbone->result.effective_stress,
              stress_step,
              [](const SymmetricTensor & stress) {
                return AbaqusCDPFormula::stressInvariants(stress).tension_weight;
              });
          const auto tensile_measure_gradient = scalarTensorGradient(
              plastic_increment, plastic_step, positivePrincipalMeasure);
          const auto compressive_measure_gradient = scalarTensorGradient(
              plastic_increment, plastic_step, negativePrincipalMeasure);

          for (std::size_t input = 0;
               input < AbaqusCDPLocalIntegrator::transition_size;
               ++input)
          {
            SymmetricTensor effective_stress_derivative = {};
            SymmetricTensor plastic_increment_derivative = {};
            for (std::size_t component = 0; component < 6; ++component)
            {
              effective_stress_derivative[component] =
                  final_backbone->derivative[input][component];
              plastic_increment_derivative[component] =
                  final_backbone->derivative[input][6 + component] -
                  (input == 6 + component ? 1.0 : 0.0);
            }

            const double weight_derivative =
                tensorDot(tension_weight_gradient, effective_stress_derivative);
            const double tensile_measure_derivative =
                tensorDot(tensile_measure_gradient, plastic_increment_derivative);
            const double compressive_measure_derivative =
                tensorDot(compressive_measure_gradient, plastic_increment_derivative);
            final_backbone->derivative[input][12] =
                (input == 12 ? 1.0 : 0.0) + weight_derivative * tensile_measure +
                tension_weight * tensile_measure_derivative;
            final_backbone->derivative[input][13] =
                (input == 13 ? 1.0 : 0.0) - weight_derivative * compressive_measure +
                (1.0 - tension_weight) * compressive_measure_derivative;
          }

          final_backbone->result = _state_integrator.aggregateBackboneHistory(
              old_state.backbone, final_backbone->result);
          if (CDPDiagnostics::current && CDPDiagnostics::current->trace)
          {
            const auto plastic_increment_principal =
                AbaqusCDPFormula::stressInvariants(plastic_increment).principal_stress;
            const auto effective_invariants =
                AbaqusCDPFormula::stressInvariants(final_backbone->result.effective_stress);
            const double raw_delta_kappa_t = tension_weight * tensile_measure;
            const double raw_delta_kappa_c = (1.0 - tension_weight) * compressive_measure;
            std::ostringstream payload;
            payload << "\"history_integration_rule\":\"global_accepted_step_end\"";
            payload << ",\"effective_stress\":";
            CDPDiagnostics::json(payload, final_backbone->result.effective_stress);
            payload << ",\"effective_principal\":";
            CDPDiagnostics::json(payload, effective_invariants.principal_stress);
            payload << ",\"plastic_increment\":";
            CDPDiagnostics::json(payload, plastic_increment);
            payload << ",\"plastic_increment_principal\":";
            CDPDiagnostics::json(payload, plastic_increment_principal);
            payload << ",\"tension_weight\":";
            CDPDiagnostics::json(payload, tension_weight);
            payload << ",\"tensile_measure\":";
            CDPDiagnostics::json(payload, tensile_measure);
            payload << ",\"compressive_measure\":";
            CDPDiagnostics::json(payload, compressive_measure);
            payload << ",\"raw_delta_kappa_t\":";
            CDPDiagnostics::json(payload, raw_delta_kappa_t);
            payload << ",\"raw_delta_kappa_c\":";
            CDPDiagnostics::json(payload, raw_delta_kappa_c);
            payload << ",\"old_kappa_t\":";
            CDPDiagnostics::json(payload, old_state.backbone.tensile_equivalent_plastic_strain);
            payload << ",\"new_kappa_t\":";
            CDPDiagnostics::json(
                payload, final_backbone->result.state.tensile_equivalent_plastic_strain);
            payload << ",\"old_kappa_c\":";
            CDPDiagnostics::json(payload, old_state.backbone.compressive_equivalent_plastic_strain);
            payload << ",\"new_kappa_c\":";
            CDPDiagnostics::json(
                payload, final_backbone->result.state.compressive_equivalent_plastic_strain);
            payload << ",\"scalable_algorithmic_tangent\":true";
            CDPDiagnostics::event("accepted_step_history", payload.str());
          }
        }
        auto assembled = _state_integrator.assembleBackboneLinearized(
            new_total_strain, time_step, old_state, *final_backbone);
        TangentMatrix tangent = {};
        for (std::size_t column = 0; column < 6; ++column)
          for (std::size_t row = 0; row < 6; ++row)
            tangent[column][row] = assembled.derivative[column][row];

        Result result{std::move(assembled.result),
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
      catch (const std::exception & error)
      {
        partition_succeeded = false;
        last_error = error.what();
        last_failed_substep = substeps;
        last_partition = substeps;
      }
    }

    partition_scope.failed();
    if (substeps == _parameters.maximum_substeps)
      break;
    substeps *= 2;
    ++cutback_count;
  }

  std::ostringstream message;
  message << "deferred-viscous linearization failed after partition=" << last_partition
          << ", failed_substep=" << last_failed_substep
          << ", attempted_partitions=" << attempted_partitions
          << ", last_error=" << last_error;
  substepError(message.str());
}

AbaqusCDPSubstepIntegrator::LinearizedResult
AbaqusCDPSubstepIntegrator::integrateAggregateHistoryReferenceLinearized(
    const SymmetricTensor & old_total_strain,
    const SymmetricTensor & new_total_strain,
    const double time_step,
    const State & old_state) const
{
  auto result = integrate(old_total_strain, new_total_strain, time_step, old_state);
  const auto reference =
      referenceTangent(old_total_strain, new_total_strain, time_step, old_state);
  return {std::move(result), reference.value};
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

  if (_parameters.aggregate_backbone_history_to_global_step &&
      _parameters.integrate_backbone_history_weights_over_substeps)
    return integrateAggregateHistoryReferenceLinearized(
        old_total_strain, new_total_strain, time_step, old_state);

  if (_parameters.defer_viscous_update_to_global_step)
    return integrateDeferredViscousLinearized(
        old_total_strain, new_total_strain, time_step, old_state);

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
