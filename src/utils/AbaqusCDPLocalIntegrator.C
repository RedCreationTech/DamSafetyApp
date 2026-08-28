#include "AbaqusCDPLocalIntegrator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace
{
using SymmetricTensor = AbaqusCDPLocalIntegrator::SymmetricTensor;

[[noreturn]] void
integrationError(const std::string & message)
{
  throw std::runtime_error("AbaqusCDPLocalIntegrator: " + message);
}

void
requireIntegratorFinite(const double value, const std::string & name)
{
  if (!std::isfinite(value))
    integrationError(name + " must be finite");
}

template <typename Vector>
double
infinityNorm(const Vector & vector)
{
  double result = 0.0;
  for (const double value : vector)
    result = std::max(result, std::abs(value));
  return result;
}

SymmetricTensor
subtract(const SymmetricTensor & left, const SymmetricTensor & right)
{
  SymmetricTensor result;
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = left[i] - right[i];
  return result;
}

SymmetricTensor
add(const SymmetricTensor & left, const SymmetricTensor & right)
{
  SymmetricTensor result;
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = left[i] + right[i];
  return result;
}

SymmetricTensor
scale(const SymmetricTensor & tensor, const double factor)
{
  SymmetricTensor result;
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = factor * tensor[i];
  return result;
}

SymmetricTensor
flowDirection(const SymmetricTensor & stress,
              const double dilation_angle_degrees,
              const double eccentricity,
              const double initial_tension_strength)
{
  auto direction = AbaqusCDPFormula::flowPotential(
                       stress, dilation_angle_degrees, eccentricity, initial_tension_strength)
                       .gradient;
  // B-006A exposes derivatives for independent symmetric components. Plastic
  // strain uses tensor shear, so the off-diagonal tensor entries are half.
  direction[3] *= 0.5;
  direction[4] *= 0.5;
  direction[5] *= 0.5;
  return direction;
}

bool
finiteTensor(const SymmetricTensor & tensor)
{
  return std::all_of(
      tensor.begin(), tensor.end(), [](const double value) { return std::isfinite(value); });
}
}

AbaqusCDPLocalIntegrator::AbaqusCDPLocalIntegrator(const CDPMaterialTable & table,
                                                   Parameters parameters)
  : _table(table),
    _parameters(parameters),
    _initial_tension_strength(
        table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0.0)
            .stress.value),
    _initial_compression_strength(
        table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::COMPRESSION, 0.0)
            .stress.value),
    _stress_scale_floor(std::max({_initial_tension_strength, _initial_compression_strength, 1.0})),
    _strain_scale(std::max(
        {table.equivalentPlasticStrain(CDPMaterialTable::Branch::TENSION).back(),
         table.equivalentPlasticStrain(CDPMaterialTable::Branch::COMPRESSION).back(),
         1.0e-6}))
{
  requireIntegratorFinite(_parameters.youngs_modulus, "Young's modulus");
  requireIntegratorFinite(_parameters.poissons_ratio, "Poisson's ratio");
  requireIntegratorFinite(_parameters.residual_tolerance, "residual tolerance");
  requireIntegratorFinite(_parameters.finite_difference_step, "finite-difference step");
  requireIntegratorFinite(_parameters.minimum_line_search, "minimum line-search factor");
  if (_parameters.youngs_modulus <= 0.0)
    integrationError("Young's modulus must be positive");
  if (_parameters.poissons_ratio <= -1.0 || _parameters.poissons_ratio >= 0.5)
    integrationError("Poisson's ratio must be in (-1, 0.5)");
  if (_parameters.residual_tolerance <= 0.0 || _parameters.finite_difference_step <= 0.0)
    integrationError("Newton tolerance and finite-difference step must be positive");
  if (_parameters.minimum_line_search <= 0.0 || _parameters.minimum_line_search > 1.0)
    integrationError("minimum line-search factor must be in (0, 1]");

  AbaqusCDPFormula::yieldCoefficients(_parameters.biaxial_to_uniaxial_compression_ratio,
                                      _parameters.tensile_meridian_ratio,
                                      _initial_compression_strength,
                                      _initial_tension_strength);
  AbaqusCDPFormula::flowPotential({0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                  _parameters.dilation_angle_degrees,
                                  _parameters.eccentricity,
                                  _initial_tension_strength);
}

AbaqusCDPLocalIntegrator::SymmetricTensor
AbaqusCDPLocalIntegrator::elasticStress(const SymmetricTensor & elastic_strain) const
{
  if (!finiteTensor(elastic_strain))
    integrationError("elastic strain contains a non-finite component");
  const double shear_modulus =
      _parameters.youngs_modulus / (2.0 * (1.0 + _parameters.poissons_ratio));
  const double lame_lambda = _parameters.youngs_modulus * _parameters.poissons_ratio /
                             ((1.0 + _parameters.poissons_ratio) *
                              (1.0 - 2.0 * _parameters.poissons_ratio));
  const double trace = elastic_strain[0] + elastic_strain[1] + elastic_strain[2];
  return {lame_lambda * trace + 2.0 * shear_modulus * elastic_strain[0],
          lame_lambda * trace + 2.0 * shear_modulus * elastic_strain[1],
          lame_lambda * trace + 2.0 * shear_modulus * elastic_strain[2],
          2.0 * shear_modulus * elastic_strain[3],
          2.0 * shear_modulus * elastic_strain[4],
          2.0 * shear_modulus * elastic_strain[5]};
}

AbaqusCDPLocalIntegrator::SymmetricTensor
AbaqusCDPLocalIntegrator::elasticStrain(const SymmetricTensor & stress) const
{
  if (!finiteTensor(stress))
    integrationError("stress contains a non-finite component");
  const double trace = stress[0] + stress[1] + stress[2];
  const double normal_factor =
      (1.0 + _parameters.poissons_ratio) / _parameters.youngs_modulus;
  const double trace_factor = _parameters.poissons_ratio / _parameters.youngs_modulus;
  const double shear_factor =
      (1.0 + _parameters.poissons_ratio) / _parameters.youngs_modulus;
  return {normal_factor * stress[0] - trace_factor * trace,
          normal_factor * stress[1] - trace_factor * trace,
          normal_factor * stress[2] - trace_factor * trace,
          shear_factor * stress[3],
          shear_factor * stress[4],
          shear_factor * stress[5]};
}

CDPMaterialTable::Response
AbaqusCDPLocalIntegrator::materialResponse(const CDPMaterialTable::Branch branch,
                                           const double equivalent_plastic_strain) const
{
  return _table.responseByEquivalentPlasticStrain(branch, equivalent_plastic_strain);
}

double
AbaqusCDPLocalIntegrator::stressScale(const SymmetricTensor & total_strain,
                                      const State & old_state) const
{
  const auto trial_stress = elasticStress(subtract(total_strain, old_state.plastic_strain));
  double result = _stress_scale_floor;
  for (const double value : trial_stress)
    result = std::max(result, std::abs(value));
  return result;
}

AbaqusCDPLocalIntegrator::Evaluation
AbaqusCDPLocalIntegrator::evaluate(const LocalVector & unknown,
                                   const SymmetricTensor & total_strain,
                                   const State & old_state,
                                   const double stress_scale) const
{
  Evaluation result;
  const auto trial_stress = elasticStress(subtract(total_strain, old_state.plastic_strain));
  for (std::size_t i = 0; i < 6; ++i)
    result.stress[i] = unknown[i] * stress_scale;
  result.plastic_multiplier = unknown[6] * _strain_scale;
  result.tensile_equivalent_plastic_strain =
      old_state.tensile_equivalent_plastic_strain + unknown[7] * _strain_scale;
  result.compressive_equivalent_plastic_strain =
      old_state.compressive_equivalent_plastic_strain + unknown[8] * _strain_scale;

  const auto direction = flowDirection(result.stress,
                                       _parameters.dilation_angle_degrees,
                                       _parameters.eccentricity,
                                       _initial_tension_strength);
  result.plastic_increment = scale(direction, result.plastic_multiplier);
  const auto elastic_correction = elasticStress(result.plastic_increment);
  for (std::size_t i = 0; i < 6; ++i)
    result.residual[i] =
        (result.stress[i] - trial_stress[i] + elastic_correction[i]) / stress_scale;

  const auto stress_invariants = AbaqusCDPFormula::stressInvariants(result.stress);
  const auto increment_principal =
      AbaqusCDPFormula::stressInvariants(result.plastic_increment).principal_stress;
  result.tensile_increment = stress_invariants.tension_weight *
                             std::max(0.0, increment_principal.back());
  result.compressive_increment = (1.0 - stress_invariants.tension_weight) *
                                 std::max(0.0, -increment_principal.front());

  const auto tension = _table.responseByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::TENSION, result.tensile_equivalent_plastic_strain);
  const auto compression = _table.responseByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::COMPRESSION, result.compressive_equivalent_plastic_strain);
  result.yield = AbaqusCDPFormula::yieldFunction(
      result.stress,
      compression.stress.value,
      tension.stress.value,
      _parameters.biaxial_to_uniaxial_compression_ratio,
      _parameters.tensile_meridian_ratio);
  result.residual[6] = result.yield / stress_scale;
  result.residual[7] =
      (result.tensile_equivalent_plastic_strain -
       old_state.tensile_equivalent_plastic_strain - result.tensile_increment) /
      _strain_scale;
  result.residual[8] =
      (result.compressive_equivalent_plastic_strain -
       old_state.compressive_equivalent_plastic_strain - result.compressive_increment) /
      _strain_scale;
  return result;
}

AbaqusCDPLocalIntegrator::LocalMatrix
AbaqusCDPLocalIntegrator::numericalJacobian(const LocalVector & unknown,
                                            const SymmetricTensor & total_strain,
                                            const State & old_state,
                                            const double stress_scale) const
{
  LocalMatrix jacobian;
  for (std::size_t column = 0; column < local_size; ++column)
  {
    const double step =
        _parameters.finite_difference_step * std::max(1.0, std::abs(unknown[column]));
    auto plus = unknown;
    auto minus = unknown;
    plus[column] += step;
    minus[column] -= step;
    const auto plus_residual = evaluate(plus, total_strain, old_state, stress_scale).residual;
    const auto minus_residual = evaluate(minus, total_strain, old_state, stress_scale).residual;
    for (std::size_t row = 0; row < local_size; ++row)
      jacobian[row][column] =
          (plus_residual[row] - minus_residual[row]) / (2.0 * step);
  }
  return jacobian;
}

AbaqusCDPLocalIntegrator::LocalVector
AbaqusCDPLocalIntegrator::solveLinearSystem(LocalMatrix matrix, LocalVector right_hand_side)
{
  for (std::size_t column = 0; column < local_size; ++column)
  {
    std::size_t pivot = column;
    for (std::size_t row = column + 1; row < local_size; ++row)
      if (std::abs(matrix[row][column]) > std::abs(matrix[pivot][column]))
        pivot = row;

    if (std::abs(matrix[pivot][column]) < 1.0e-14)
      integrationError("singular numerical local Jacobian");
    if (pivot != column)
    {
      std::swap(matrix[pivot], matrix[column]);
      std::swap(right_hand_side[pivot], right_hand_side[column]);
    }

    const double diagonal = matrix[column][column];
    for (std::size_t entry = column; entry < local_size; ++entry)
      matrix[column][entry] /= diagonal;
    right_hand_side[column] /= diagonal;

    for (std::size_t row = 0; row < local_size; ++row)
    {
      if (row == column)
        continue;
      const double factor = matrix[row][column];
      for (std::size_t entry = column; entry < local_size; ++entry)
        matrix[row][entry] -= factor * matrix[column][entry];
      right_hand_side[row] -= factor * right_hand_side[column];
    }
  }
  return right_hand_side;
}

AbaqusCDPLocalIntegrator::Result
AbaqusCDPLocalIntegrator::integrate(const SymmetricTensor & total_strain,
                                    const State & old_state) const
{
  if (!finiteTensor(total_strain) || !finiteTensor(old_state.plastic_strain))
    integrationError("total or plastic strain contains a non-finite component");
  requireIntegratorFinite(old_state.tensile_equivalent_plastic_strain,
                          "old tensile equivalent plastic strain");
  requireIntegratorFinite(old_state.compressive_equivalent_plastic_strain,
                          "old compressive equivalent plastic strain");
  if (old_state.tensile_equivalent_plastic_strain < 0.0 ||
      old_state.compressive_equivalent_plastic_strain < 0.0)
    integrationError("old equivalent plastic strains must be nonnegative");

  const auto trial_stress = elasticStress(subtract(total_strain, old_state.plastic_strain));
  const auto old_tension = _table.responseByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::TENSION, old_state.tensile_equivalent_plastic_strain);
  const auto old_compression = _table.responseByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::COMPRESSION, old_state.compressive_equivalent_plastic_strain);
  const double stress_scale = stressScale(total_strain, old_state);

  const double trial_yield = AbaqusCDPFormula::yieldFunction(
      trial_stress,
      old_compression.stress.value,
      old_tension.stress.value,
      _parameters.biaxial_to_uniaxial_compression_ratio,
      _parameters.tensile_meridian_ratio);
  if (trial_yield <= _parameters.residual_tolerance * stress_scale)
    return {trial_stress,
            old_state,
            ActiveBranch::ELASTIC,
            false,
            0,
            0.0,
            trial_yield,
            trial_yield,
            0.0,
            old_tension.damage.value,
            old_compression.damage.value};

  LocalVector unknown = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double initial_multiplier = std::max(trial_yield / _parameters.youngs_modulus, 1.0e-14);
  const auto initial_direction = flowDirection(trial_stress,
                                               _parameters.dilation_angle_degrees,
                                               _parameters.eccentricity,
                                               _initial_tension_strength);
  const auto initial_increment = scale(initial_direction, initial_multiplier);
  const auto initial_stress = subtract(trial_stress, elasticStress(initial_increment));
  for (std::size_t i = 0; i < 6; ++i)
    unknown[i] = initial_stress[i] / stress_scale;
  unknown[6] = initial_multiplier / _strain_scale;
  const auto initial_stress_invariants = AbaqusCDPFormula::stressInvariants(initial_stress);
  const auto initial_increment_principal =
      AbaqusCDPFormula::stressInvariants(initial_increment).principal_stress;
  unknown[7] = initial_stress_invariants.tension_weight *
               std::max(0.0, initial_increment_principal.back()) / _strain_scale;
  unknown[8] = (1.0 - initial_stress_invariants.tension_weight) *
               std::max(0.0, -initial_increment_principal.front()) / _strain_scale;

  Evaluation current = evaluate(unknown, total_strain, old_state, stress_scale);
  double residual_norm = infinityNorm(current.residual);
  const auto currentBranch = [&](const Evaluation & evaluation) {
    const double tolerance = 10.0 * std::numeric_limits<double>::epsilon() * _strain_scale;
    if (evaluation.tensile_increment <= tolerance)
      return ActiveBranch::COMPRESSION;
    if (evaluation.compressive_increment <= tolerance)
      return ActiveBranch::TENSION;
    return ActiveBranch::MIXED;
  };
  unsigned int iterations = 0;
  for (; iterations < _parameters.maximum_iterations &&
         residual_norm > _parameters.residual_tolerance;
       ++iterations)
  {
    const auto jacobian = numericalJacobian(unknown, total_strain, old_state, stress_scale);

    LocalVector right_hand_side;
    for (std::size_t i = 0; i < local_size; ++i)
      right_hand_side[i] = -current.residual[i];
    const auto increment = solveLinearSystem(jacobian, right_hand_side);

    bool accepted = false;
    for (double line_search = 1.0; line_search >= _parameters.minimum_line_search;
         line_search *= 0.5)
    {
      auto candidate_unknown = unknown;
      for (std::size_t i = 0; i < local_size; ++i)
        candidate_unknown[i] += line_search * increment[i];
      if (candidate_unknown[6] < 0.0 || candidate_unknown[7] < 0.0 ||
          candidate_unknown[8] < 0.0)
        continue;

      auto candidate = evaluate(candidate_unknown, total_strain, old_state, stress_scale);
      const double candidate_norm = infinityNorm(candidate.residual);
      if (std::isfinite(candidate_norm) && candidate_norm < residual_norm)
      {
        unknown = candidate_unknown;
        current = candidate;
        residual_norm = candidate_norm;
        accepted = true;
        break;
      }
    }
    if (!accepted)
    {
      std::ostringstream message;
      message << std::setprecision(16) << "line search failed at iteration " << iterations + 1
              << ", residual=" << residual_norm
              << ", plastic_multiplier=" << current.plastic_multiplier
              << ", kappa_t=" << current.tensile_equivalent_plastic_strain
              << ", kappa_c=" << current.compressive_equivalent_plastic_strain
              << ", branch=" << branchName(currentBranch(current));
      integrationError(message.str());
    }
  }

  if (residual_norm > _parameters.residual_tolerance)
  {
    std::ostringstream message;
    message << std::setprecision(16) << "local Newton failed after " << iterations
            << " iterations, residual=" << residual_norm
            << ", plastic_multiplier=" << current.plastic_multiplier
            << ", kappa_t=" << current.tensile_equivalent_plastic_strain
            << ", kappa_c=" << current.compressive_equivalent_plastic_strain
            << ", branch=" << branchName(currentBranch(current));
    integrationError(message.str());
  }

  State new_state = old_state;
  new_state.plastic_strain = add(old_state.plastic_strain, current.plastic_increment);
  new_state.tensile_equivalent_plastic_strain = current.tensile_equivalent_plastic_strain;
  new_state.compressive_equivalent_plastic_strain =
      current.compressive_equivalent_plastic_strain;
  const auto tension = _table.responseByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::TENSION, new_state.tensile_equivalent_plastic_strain);
  const auto compression = _table.responseByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::COMPRESSION, new_state.compressive_equivalent_plastic_strain);

  return {current.stress,
          new_state,
          currentBranch(current),
          true,
          iterations,
          residual_norm,
          trial_yield,
          current.yield,
          current.plastic_multiplier,
          tension.damage.value,
          compression.damage.value};
}

AbaqusCDPLocalIntegrator::LinearizedResult
AbaqusCDPLocalIntegrator::integrateLinearized(const SymmetricTensor & total_strain,
                                              const State & old_state) const
{
  LinearizedResult linearized{integrate(total_strain, old_state), {}};

  // The first six columns of the elastic stiffness use the same physical
  // tensor-shear convention as the constitutive arrays.
  std::array<SymmetricTensor, 6> elastic_columns;
  for (std::size_t column = 0; column < 6; ++column)
  {
    SymmetricTensor unit = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    unit[column] = 1.0;
    elastic_columns[column] = elasticStress(unit);
  }

  if (!linearized.result.plastic)
  {
    for (std::size_t column = 0; column < 6; ++column)
      for (std::size_t row = 0; row < 6; ++row)
      {
        linearized.derivative[column][row] = elastic_columns[column][row];
        linearized.derivative[6 + column][row] = -elastic_columns[column][row];
        linearized.derivative[6 + column][6 + row] = column == row ? 1.0 : 0.0;
      }
    linearized.derivative[12][12] = 1.0;
    linearized.derivative[13][13] = 1.0;
    return linearized;
  }

  const double stress_scale = stressScale(total_strain, old_state);
  LocalVector unknown = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (std::size_t i = 0; i < 6; ++i)
    unknown[i] = linearized.result.effective_stress[i] / stress_scale;
  unknown[6] = linearized.result.plastic_multiplier / _strain_scale;
  unknown[7] = (linearized.result.state.tensile_equivalent_plastic_strain -
                old_state.tensile_equivalent_plastic_strain) /
               _strain_scale;
  unknown[8] = (linearized.result.state.compressive_equivalent_plastic_strain -
                old_state.compressive_equivalent_plastic_strain) /
               _strain_scale;
  const auto jacobian = numericalJacobian(unknown, total_strain, old_state, stress_scale);

  const auto tension = materialResponse(
      CDPMaterialTable::Branch::TENSION,
      linearized.result.state.tensile_equivalent_plastic_strain);
  const auto compression = materialResponse(
      CDPMaterialTable::Branch::COMPRESSION,
      linearized.result.state.compressive_equivalent_plastic_strain);
  const auto yield_with_strengths = [&](const double compression_strength,
                                        const double tension_strength) {
    return AbaqusCDPFormula::yieldFunction(
        linearized.result.effective_stress,
        compression_strength,
        tension_strength,
        _parameters.biaxial_to_uniaxial_compression_ratio,
        _parameters.tensile_meridian_ratio);
  };
  const double tension_step = 1.0e-7 * std::max(1.0, std::abs(tension.stress.value));
  const double compression_step = 1.0e-7 * std::max(1.0, std::abs(compression.stress.value));
  const double yield_tension_derivative =
      (yield_with_strengths(compression.stress.value, tension.stress.value + tension_step) -
       yield_with_strengths(compression.stress.value, tension.stress.value - tension_step)) /
      (2.0 * tension_step);
  const double yield_compression_derivative =
      (yield_with_strengths(compression.stress.value + compression_step, tension.stress.value) -
       yield_with_strengths(compression.stress.value - compression_step, tension.stress.value)) /
      (2.0 * compression_step);

  for (std::size_t input = 0; input < transition_size; ++input)
  {
    LocalVector direct = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (input < 6)
      for (std::size_t row = 0; row < 6; ++row)
        direct[row] = -elastic_columns[input][row] / stress_scale;
    else if (input < 12)
      for (std::size_t row = 0; row < 6; ++row)
        direct[row] = elastic_columns[input - 6][row] / stress_scale;
    else if (input == 12)
      direct[6] =
          yield_tension_derivative * tension.stress.right_derivative / stress_scale;
    else
      direct[6] =
          yield_compression_derivative * compression.stress.right_derivative / stress_scale;

    for (double & value : direct)
      value = -value;
    const auto unknown_derivative = solveLinearSystem(jacobian, direct);

    SymmetricTensor stress_derivative;
    for (std::size_t row = 0; row < 6; ++row)
    {
      stress_derivative[row] = stress_scale * unknown_derivative[row];
      linearized.derivative[input][row] = stress_derivative[row];
    }
    const auto elastic_strain_derivative = elasticStrain(stress_derivative);
    for (std::size_t row = 0; row < 6; ++row)
      linearized.derivative[input][6 + row] =
          (input == row ? 1.0 : 0.0) - elastic_strain_derivative[row];
    linearized.derivative[input][12] =
        (input == 12 ? 1.0 : 0.0) + _strain_scale * unknown_derivative[7];
    linearized.derivative[input][13] =
        (input == 13 ? 1.0 : 0.0) + _strain_scale * unknown_derivative[8];
  }
  return linearized;
}

std::string
AbaqusCDPLocalIntegrator::branchName(const ActiveBranch branch)
{
  switch (branch)
  {
    case ActiveBranch::ELASTIC:
      return "elastic";
    case ActiveBranch::TENSION:
      return "tension";
    case ActiveBranch::COMPRESSION:
      return "compression";
    case ActiveBranch::MIXED:
      return "mixed";
  }
  return "unknown";
}
