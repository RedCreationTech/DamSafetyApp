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
constexpr std::size_t local_size = 9;
using LocalVector = std::array<double, local_size>;
using LocalMatrix = std::array<LocalVector, local_size>;
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

double
infinityNorm(const LocalVector & vector)
{
  double result = 0.0;
  for (const double value : vector)
    result = std::max(result, std::abs(value));
  return result;
}

LocalVector
solveLinearSystem(LocalMatrix matrix, LocalVector right_hand_side)
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
  double stress_scale = _stress_scale_floor;
  for (const double value : trial_stress)
    stress_scale = std::max(stress_scale, std::abs(value));

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

  struct Evaluation
  {
    LocalVector residual;
    SymmetricTensor stress;
    SymmetricTensor plastic_increment;
    double plastic_multiplier;
    double tensile_equivalent_plastic_strain;
    double compressive_equivalent_plastic_strain;
    double tensile_increment;
    double compressive_increment;
    double yield;
  };

  const auto evaluate = [&](const LocalVector & unknown) {
    Evaluation result;
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
  };

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

  Evaluation current = evaluate(unknown);
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
    LocalMatrix jacobian;
    for (std::size_t column = 0; column < local_size; ++column)
    {
      const double step =
          _parameters.finite_difference_step * std::max(1.0, std::abs(unknown[column]));
      auto plus = unknown;
      auto minus = unknown;
      plus[column] += step;
      minus[column] -= step;
      const auto plus_residual = evaluate(plus).residual;
      const auto minus_residual = evaluate(minus).residual;
      for (std::size_t row = 0; row < local_size; ++row)
        jacobian[row][column] =
            (plus_residual[row] - minus_residual[row]) / (2.0 * step);
    }

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

      auto candidate = evaluate(candidate_unknown);
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
