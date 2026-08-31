#include "CDPDiagnostics.h"
#include "AbaqusCDPLocalIntegrator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace
{
using SymmetricTensor = AbaqusCDPLocalIntegrator::SymmetricTensor;
constexpr std::size_t ad_size = AbaqusCDPLocalIntegrator::local_size;

struct Dual9
{
  double value = 0.0;
  std::array<double, ad_size> derivative = {};

  Dual9() = default;
  Dual9(const double input) : value(input) {}

  static Dual9 variable(const double input, const std::size_t index, const double scale = 1.0)
  {
    Dual9 result(input);
    result.derivative[index] = scale;
    return result;
  }
};

Dual9
operator+(const Dual9 & left, const Dual9 & right)
{
  Dual9 result(left.value + right.value);
  for (std::size_t i = 0; i < ad_size; ++i)
    result.derivative[i] = left.derivative[i] + right.derivative[i];
  return result;
}

Dual9
operator-(const Dual9 & left, const Dual9 & right)
{
  Dual9 result(left.value - right.value);
  for (std::size_t i = 0; i < ad_size; ++i)
    result.derivative[i] = left.derivative[i] - right.derivative[i];
  return result;
}

Dual9
operator-(const Dual9 & input)
{
  Dual9 result(-input.value);
  for (std::size_t i = 0; i < ad_size; ++i)
    result.derivative[i] = -input.derivative[i];
  return result;
}

Dual9
operator*(const Dual9 & left, const Dual9 & right)
{
  Dual9 result(left.value * right.value);
  for (std::size_t i = 0; i < ad_size; ++i)
    result.derivative[i] =
        left.derivative[i] * right.value + left.value * right.derivative[i];
  return result;
}

Dual9
operator/(const Dual9 & numerator, const Dual9 & denominator)
{
  Dual9 result(numerator.value / denominator.value);
  const double denominator_squared = denominator.value * denominator.value;
  for (std::size_t i = 0; i < ad_size; ++i)
    result.derivative[i] =
        (numerator.derivative[i] * denominator.value -
         numerator.value * denominator.derivative[i]) /
        denominator_squared;
  return result;
}

Dual9
dualSqrt(const Dual9 & input)
{
  if (input.value <= 0.0)
    return Dual9(0.0);
  Dual9 result(std::sqrt(input.value));
  for (std::size_t i = 0; i < ad_size; ++i)
    result.derivative[i] = input.derivative[i] / (2.0 * result.value);
  return result;
}

Dual9
dualAcos(const Dual9 & input)
{
  if (input.value <= -1.0)
    return Dual9(std::acos(-1.0));
  if (input.value >= 1.0)
    return Dual9(0.0);
  Dual9 result(std::acos(input.value));
  const double factor = -1.0 / std::sqrt(1.0 - input.value * input.value);
  for (std::size_t i = 0; i < ad_size; ++i)
    result.derivative[i] = factor * input.derivative[i];
  return result;
}

Dual9
dualCos(const Dual9 & input)
{
  Dual9 result(std::cos(input.value));
  const double factor = -std::sin(input.value);
  for (std::size_t i = 0; i < ad_size; ++i)
    result.derivative[i] = factor * input.derivative[i];
  return result;
}

Dual9
dualPositivePart(const Dual9 & input)
{
  return input.value > 0.0 ? input : Dual9(0.0);
}

Dual9
dualAbsoluteValue(const Dual9 & input)
{
  if (input.value > 0.0)
    return input;
  if (input.value < 0.0)
    return -input;
  return Dual9(0.0);
}

using DualTensor = std::array<Dual9, 6>;

std::array<Dual9, 3>
dualPrincipalStress(const DualTensor & stress)
{
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::AD_SPECTRUM);
  constexpr double pi = 3.141592653589793238462643383279502884;
  const Dual9 & xx = stress[0];
  const Dual9 & yy = stress[1];
  const Dual9 & zz = stress[2];
  const Dual9 & xy = stress[3];
  const Dual9 & yz = stress[4];
  const Dual9 & xz = stress[5];
  const Dual9 off_diagonal_norm = xy * xy + yz * yz + xz * xz;

  std::array<Dual9, 3> result;
  if (off_diagonal_norm.value == 0.0)
    result = {xx, yy, zz};
  else
  {
    const Dual9 center = (xx + yy + zz) / 3.0;
    const Dual9 centered_norm = (xx - center) * (xx - center) +
                                (yy - center) * (yy - center) +
                                (zz - center) * (zz - center) + 2.0 * off_diagonal_norm;
    const Dual9 scale = dualSqrt(centered_norm / 6.0);
    if (scale.value == 0.0)
      result = {center, center, center};
    else
    {
      const Dual9 bxx = (xx - center) / scale;
      const Dual9 byy = (yy - center) / scale;
      const Dual9 bzz = (zz - center) / scale;
      const Dual9 bxy = xy / scale;
      const Dual9 byz = yz / scale;
      const Dual9 bxz = xz / scale;
      const Dual9 determinant = bxx * byy * bzz + 2.0 * bxy * byz * bxz -
                                bxx * byz * byz - byy * bxz * bxz - bzz * bxy * bxy;
      const Dual9 angle = dualAcos(determinant / 2.0) / 3.0;
      const Dual9 largest = center + 2.0 * scale * dualCos(angle);
      const Dual9 smallest = center + 2.0 * scale * dualCos(angle + 2.0 * pi / 3.0);
      result = {smallest, 3.0 * center - largest - smallest, largest};
    }
  }

  std::sort(result.begin(), result.end(), [](const Dual9 & left, const Dual9 & right) {
    return left.value < right.value;
  });
  return result;
}

struct DualStressInvariants
{
  Dual9 pressure;
  Dual9 mises;
  std::array<Dual9, 3> principal_stress;
  Dual9 maximum_principal_stress;
  Dual9 tension_weight;
};

DualStressInvariants
dualStressInvariants(const DualTensor & stress)
{
  const Dual9 trace = stress[0] + stress[1] + stress[2];
  const Dual9 mean = trace / 3.0;
  const Dual9 sxx = stress[0] - mean;
  const Dual9 syy = stress[1] - mean;
  const Dual9 szz = stress[2] - mean;
  const Dual9 j2 = 0.5 * (sxx * sxx + syy * syy + szz * szz +
                          2.0 * (stress[3] * stress[3] + stress[4] * stress[4] +
                                 stress[5] * stress[5]));
  const auto principal = dualPrincipalStress(stress);
  Dual9 positive_sum;
  Dual9 absolute_sum;
  for (const auto & value : principal)
  {
    positive_sum = positive_sum + dualPositivePart(value);
    absolute_sum = absolute_sum + dualAbsoluteValue(value);
  }
  return {-mean,
          dualSqrt(dualPositivePart(3.0 * j2)),
          principal,
          principal.back(),
          absolute_sum.value == 0.0 ? Dual9(0.0) : positive_sum / absolute_sum};
}

DualTensor
dualFlowDirection(const DualTensor & stress,
                  const double dilation_angle_degrees,
                  const double eccentricity,
                  const double initial_tension_strength)
{
  constexpr double pi = 3.141592653589793238462643383279502884;
  const auto invariants = dualStressInvariants(stress);
  const double tangent = std::tan(dilation_angle_degrees * pi / 180.0);
  const double asymptote = eccentricity * initial_tension_strength * tangent;
  const Dual9 root = dualSqrt(asymptote * asymptote + invariants.mises * invariants.mises);
  const Dual9 mean = (stress[0] + stress[1] + stress[2]) / 3.0;

  DualTensor direction;
  direction[0] = 1.5 * (stress[0] - mean) / root + tangent / 3.0;
  direction[1] = 1.5 * (stress[1] - mean) / root + tangent / 3.0;
  direction[2] = 1.5 * (stress[2] - mean) / root + tangent / 3.0;
  // Plastic strain stores physical tensor shear components, half of the
  // ordinary derivative with respect to an independent symmetric shear.
  direction[3] = 1.5 * stress[3] / root;
  direction[4] = 1.5 * stress[4] / root;
  direction[5] = 1.5 * stress[5] / root;
  return direction;
}

Dual9
dualYieldFunction(const DualTensor & stress,
                  const Dual9 & compression_strength,
                  const Dual9 & tension_strength,
                  const double biaxial_to_uniaxial_compression_ratio,
                  const double tensile_meridian_ratio)
{
  const auto invariants = dualStressInvariants(stress);
  const double alpha = (biaxial_to_uniaxial_compression_ratio - 1.0) /
                       (2.0 * biaxial_to_uniaxial_compression_ratio - 1.0);
  const Dual9 beta = compression_strength / tension_strength * (1.0 - alpha) -
                     (1.0 + alpha);
  const double gamma = 3.0 * (1.0 - tensile_meridian_ratio) /
                       (2.0 * tensile_meridian_ratio - 1.0);
  const Dual9 maximum = invariants.maximum_principal_stress;
  return (invariants.mises - 3.0 * alpha * invariants.pressure +
          beta * dualPositivePart(maximum) - gamma * dualPositivePart(-maximum)) /
             (1.0 - alpha) -
         compression_strength;
}

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
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::RESIDUAL);
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

  const auto tension = _table.effectiveStrengthByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::TENSION, result.tensile_equivalent_plastic_strain);
  const auto compression = _table.effectiveStrengthByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::COMPRESSION, result.compressive_equivalent_plastic_strain);
  result.yield = AbaqusCDPFormula::yieldFunction(
      result.stress,
      compression.value,
      tension.value,
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
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::FD_JACOBIAN);
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

AbaqusCDPLocalIntegrator::LocalMatrix
AbaqusCDPLocalIntegrator::automaticDifferentiationJacobian(const LocalVector & unknown,
                                                           const SymmetricTensor & total_strain,
                                                           const State & old_state,
                                                           const double stress_scale) const
{
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::AD_JACOBIAN);
  const auto trial_stress = elasticStress(subtract(total_strain, old_state.plastic_strain));
  DualTensor stress;
  for (std::size_t i = 0; i < 6; ++i)
    stress[i] = Dual9::variable(unknown[i] * stress_scale, i, stress_scale);
  const Dual9 plastic_multiplier = Dual9::variable(unknown[6] * _strain_scale, 6, _strain_scale);
  const Dual9 tensile_equivalent_plastic_strain =
      Dual9::variable(old_state.tensile_equivalent_plastic_strain +
                          unknown[7] * _strain_scale,
                      7,
                      _strain_scale);
  const Dual9 compressive_equivalent_plastic_strain =
      Dual9::variable(old_state.compressive_equivalent_plastic_strain +
                          unknown[8] * _strain_scale,
                      8,
                      _strain_scale);

  const auto direction = dualFlowDirection(stress,
                                           _parameters.dilation_angle_degrees,
                                           _parameters.eccentricity,
                                           _initial_tension_strength);
  DualTensor plastic_increment;
  for (std::size_t i = 0; i < 6; ++i)
    plastic_increment[i] = direction[i] * plastic_multiplier;

  const double shear_modulus =
      _parameters.youngs_modulus / (2.0 * (1.0 + _parameters.poissons_ratio));
  const double lame_lambda = _parameters.youngs_modulus * _parameters.poissons_ratio /
                             ((1.0 + _parameters.poissons_ratio) *
                              (1.0 - 2.0 * _parameters.poissons_ratio));
  const Dual9 plastic_trace =
      plastic_increment[0] + plastic_increment[1] + plastic_increment[2];
  DualTensor elastic_correction;
  elastic_correction[0] =
      lame_lambda * plastic_trace + 2.0 * shear_modulus * plastic_increment[0];
  elastic_correction[1] =
      lame_lambda * plastic_trace + 2.0 * shear_modulus * plastic_increment[1];
  elastic_correction[2] =
      lame_lambda * plastic_trace + 2.0 * shear_modulus * plastic_increment[2];
  elastic_correction[3] = 2.0 * shear_modulus * plastic_increment[3];
  elastic_correction[4] = 2.0 * shear_modulus * plastic_increment[4];
  elastic_correction[5] = 2.0 * shear_modulus * plastic_increment[5];

  std::array<Dual9, local_size> residual;
  for (std::size_t i = 0; i < 6; ++i)
    residual[i] = (stress[i] - trial_stress[i] + elastic_correction[i]) / stress_scale;

  const auto stress_invariants = dualStressInvariants(stress);
  const auto increment_principal = dualPrincipalStress(plastic_increment);
  const Dual9 tensile_increment =
      stress_invariants.tension_weight * dualPositivePart(increment_principal.back());
  const Dual9 compressive_increment =
      (1.0 - stress_invariants.tension_weight) *
      dualPositivePart(-increment_principal.front());

  const auto tension = _table.effectiveStrengthByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::TENSION, tensile_equivalent_plastic_strain.value);
  const auto compression = _table.effectiveStrengthByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::COMPRESSION, compressive_equivalent_plastic_strain.value);
  const auto lift_table_sample = [](const Dual9 & abscissa,
                                    const CDPMaterialTable::Sample & sample) {
    Dual9 result(sample.value);
    for (std::size_t i = 0; i < ad_size; ++i)
      result.derivative[i] = sample.right_derivative * abscissa.derivative[i];
    return result;
  };
  const Dual9 tension_strength = lift_table_sample(tensile_equivalent_plastic_strain,
                                                   tension);
  const Dual9 compression_strength = lift_table_sample(compressive_equivalent_plastic_strain,
                                                       compression);
  residual[6] = dualYieldFunction(stress,
                                  compression_strength,
                                  tension_strength,
                                  _parameters.biaxial_to_uniaxial_compression_ratio,
                                  _parameters.tensile_meridian_ratio) /
                stress_scale;
  residual[7] = (tensile_equivalent_plastic_strain -
                 old_state.tensile_equivalent_plastic_strain - tensile_increment) /
                _strain_scale;
  residual[8] = (compressive_equivalent_plastic_strain -
                 old_state.compressive_equivalent_plastic_strain - compressive_increment) /
                _strain_scale;

  LocalMatrix jacobian;
  for (std::size_t row = 0; row < local_size; ++row)
    for (std::size_t column = 0; column < local_size; ++column)
      jacobian[row][column] = residual[row].derivative[column];
  return jacobian;
}

AbaqusCDPLocalIntegrator::LocalMatrix
AbaqusCDPLocalIntegrator::localJacobian(const LocalVector & unknown,
                                        const SymmetricTensor & total_strain,
                                        const State & old_state,
                                        const double stress_scale) const
{
  if (_parameters.use_automatic_differentiation_jacobian)
    return automaticDifferentiationJacobian(unknown, total_strain, old_state, stress_scale);
  return numericalJacobian(unknown, total_strain, old_state, stress_scale);
}

AbaqusCDPLocalIntegrator::LocalFactorization
AbaqusCDPLocalIntegrator::factorLinearSystem(LocalMatrix matrix)
{
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::FACTOR);
  LocalFactorization factorization{std::move(matrix), {}};
  for (std::size_t column = 0; column < local_size; ++column)
  {
    std::size_t pivot = column;
    for (std::size_t row = column + 1; row < local_size; ++row)
      if (std::abs(factorization.factors[row][column]) >
          std::abs(factorization.factors[pivot][column]))
        pivot = row;

    factorization.pivots[column] = pivot;
    if (std::abs(factorization.factors[pivot][column]) < 1.0e-14)
      integrationError("singular local Jacobian");
    if (pivot != column)
      std::swap(factorization.factors[pivot], factorization.factors[column]);

    const double diagonal = factorization.factors[column][column];
    for (std::size_t row = column + 1; row < local_size; ++row)
    {
      factorization.factors[row][column] /= diagonal;
      const double multiplier = factorization.factors[row][column];
      for (std::size_t entry = column + 1; entry < local_size; ++entry)
        factorization.factors[row][entry] -=
            multiplier * factorization.factors[column][entry];
    }
  }

  return factorization;
}

AbaqusCDPLocalIntegrator::LocalVector
AbaqusCDPLocalIntegrator::solveLinearSystem(const LocalFactorization & factorization,
                                            LocalVector right_hand_side)
{
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::BACKSOLVE);
  for (std::size_t column = 0; column < local_size; ++column)
    if (factorization.pivots[column] != column)
      std::swap(right_hand_side[factorization.pivots[column]], right_hand_side[column]);

  for (std::size_t row = 0; row < local_size; ++row)
    for (std::size_t column = 0; column < row; ++column)
      right_hand_side[row] -=
          factorization.factors[row][column] * right_hand_side[column];

  for (std::size_t offset = 0; offset < local_size; ++offset)
  {
    const std::size_t row = local_size - 1 - offset;
    for (std::size_t column = row + 1; column < local_size; ++column)
      right_hand_side[row] -=
          factorization.factors[row][column] * right_hand_side[column];
    right_hand_side[row] /= factorization.factors[row][row];
  }
  return right_hand_side;
}

AbaqusCDPLocalIntegrator::Result
AbaqusCDPLocalIntegrator::integrate(const SymmetricTensor & total_strain,
                                    const State & old_state) const
{
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::LOCAL);
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
  const auto old_tension = _table.effectiveStrengthByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::TENSION, old_state.tensile_equivalent_plastic_strain);
  const auto old_compression = _table.effectiveStrengthByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::COMPRESSION, old_state.compressive_equivalent_plastic_strain);
  const double stress_scale = stressScale(total_strain, old_state);

  const double trial_yield = AbaqusCDPFormula::yieldFunction(
      trial_stress,
      old_compression.value,
      old_tension.value,
      _parameters.biaxial_to_uniaxial_compression_ratio,
      _parameters.tensile_meridian_ratio);
  if (trial_yield <= _parameters.residual_tolerance * stress_scale)
    return {trial_stress,
            old_state,
            ActiveBranch::ELASTIC,
            false,
            0,
            0,
            0,
            0,
            0,
            0,
            0.0,
            trial_yield,
            trial_yield,
            0.0,
            materialResponse(CDPMaterialTable::Branch::TENSION,
                             old_state.tensile_equivalent_plastic_strain).damage.value,
            materialResponse(CDPMaterialTable::Branch::COMPRESSION,
                             old_state.compressive_equivalent_plastic_strain).damage.value};

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
  const auto diagnose_failure = [&]() {
    if (!CDPDiagnostics::sampleFailure()) return;
    std::ostringstream payload;
    try
    {
      {
        CDPDiagnostics::Binding exclude_diagnostic_cost(nullptr);
        payload<<"\"target\":";CDPDiagnostics::json(payload,total_strain);
        payload<<",\"old_plastic\":";CDPDiagnostics::json(payload,old_state.plastic_strain);
        payload<<",\"old_kappa_t\":";CDPDiagnostics::json(payload,old_state.tensile_equivalent_plastic_strain);
        payload<<",\"old_kappa_c\":";CDPDiagnostics::json(payload,old_state.compressive_equivalent_plastic_strain);
        payload<<",\"unknown\":";CDPDiagnostics::json(payload,unknown);
        payload<<",\"residual\":";CDPDiagnostics::json(payload,current.residual);
        payload<<",\"stress_scale\":";CDPDiagnostics::json(payload,stress_scale);
        payload<<",\"strain_scale\":";CDPDiagnostics::json(payload,_strain_scale);
        const auto ad=automaticDifferentiationJacobian(unknown,total_strain,old_state,stress_scale);
        payload<<",\"ad\":";CDPDiagnostics::json(payload,ad);
        payload<<",\"difference_sweep\":[";
        bool first=true;
        for (double h : {1e-5,1e-6,1e-7,1e-8})
        {
          LocalMatrix centered={}, forward={};
          for (std::size_t column=0;column<local_size;++column)
          {
            auto plus=unknown,minus=unknown;plus[column]+=h;minus[column]-=h;
            const auto rp=evaluate(plus,total_strain,old_state,stress_scale).residual;
            const auto rm=evaluate(minus,total_strain,old_state,stress_scale).residual;
            for (std::size_t row=0;row<local_size;++row)
            { centered[row][column]=(rp[row]-rm[row])/(2*h);forward[row][column]=(rp[row]-current.residual[row])/h; }
          }
          if(!first)payload<<',';first=false;
          payload<<"{\"h\":"<<h<<",\"centered\":";CDPDiagnostics::json(payload,centered);
          payload<<",\"forward\":";CDPDiagnostics::json(payload,forward);payload<<'}';
        }
        payload<<']';
      }
      CDPDiagnostics::event("local_failure_jacobian",payload.str());
    }
    catch (const std::exception &) { CDPDiagnostics::event("local_failure_jacobian","\"diagnostic_error\":true"); }
  };
  unsigned int iterations = 0;
  unsigned int jacobian_fallbacks = 0;
  unsigned int automatic_jacobian_evaluations = 0;
  unsigned int finite_difference_jacobian_evaluations = 0;
  unsigned int local_factorizations = 0;
  unsigned int local_backsolves = 0;
  for (; iterations < _parameters.maximum_iterations &&
         residual_norm > _parameters.residual_tolerance;
       ++iterations)
  {
    CDPDiagnostics::Scope newton_scope(CDPDiagnostics::NEWTON_STEP);
    LocalVector right_hand_side;
    for (std::size_t i = 0; i < local_size; ++i)
      right_hand_side[i] = -current.residual[i];
    const auto try_increment = [&](const LocalVector & increment) {
      for (double line_search = 1.0; line_search >= _parameters.minimum_line_search;
           line_search *= 0.5)
      {
        auto candidate_unknown = unknown;
        for (std::size_t i = 0; i < local_size; ++i)
          candidate_unknown[i] += line_search * increment[i];
        if (_parameters.use_bound_feasible_line_search)
        {
          if (!std::all_of(candidate_unknown.begin(), candidate_unknown.end(),
                           [](const double value) { return std::isfinite(value); }))
            continue;
          // Project trial multiplier/history increments, never the accepted history.
          // Acceptance below still requires descent in the complete original residual.
          for (std::size_t i = 6; i < local_size; ++i)
            candidate_unknown[i] = std::max(0.0, candidate_unknown[i]);
        }
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
          return true;
        }
      }
      return false;
    };

    bool accepted = false;
    if (_parameters.use_automatic_differentiation_jacobian)
    {
      try
      {
        const auto automatic = automaticDifferentiationJacobian(
            unknown, total_strain, old_state, stress_scale);
        ++automatic_jacobian_evaluations;
        const auto factorization = factorLinearSystem(automatic);
        ++local_factorizations;
        accepted = try_increment(solveLinearSystem(factorization, right_hand_side));
        ++local_backsolves;
      }
      catch (const std::runtime_error &)
      {
        accepted = false;
      }
      if (!accepted)
      {
        ++jacobian_fallbacks;
        const auto reference = numericalJacobian(unknown, total_strain, old_state, stress_scale);
        ++finite_difference_jacobian_evaluations;
        const auto factorization = factorLinearSystem(reference);
        ++local_factorizations;
        accepted = try_increment(solveLinearSystem(factorization, right_hand_side));
        ++local_backsolves;
      }
    }
    else
    {
      const auto reference = numericalJacobian(unknown, total_strain, old_state, stress_scale);
      ++finite_difference_jacobian_evaluations;
      const auto factorization = factorLinearSystem(reference);
      ++local_factorizations;
      accepted = try_increment(solveLinearSystem(factorization, right_hand_side));
      ++local_backsolves;
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
      diagnose_failure();
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
    diagnose_failure();
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
          jacobian_fallbacks,
          automatic_jacobian_evaluations,
          finite_difference_jacobian_evaluations,
          local_factorizations,
          local_backsolves,
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
  CDPDiagnostics::Scope diagnostic_scope(CDPDiagnostics::LOCAL_LINEARIZED);
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
  const auto jacobian = localJacobian(unknown, total_strain, old_state, stress_scale);
  if (_parameters.use_automatic_differentiation_jacobian)
    ++linearized.result.automatic_jacobian_evaluations;
  else
    ++linearized.result.finite_difference_jacobian_evaluations;
  const auto factorization = factorLinearSystem(jacobian);
  ++linearized.result.local_factorizations;

  const auto tension = _table.effectiveStrengthByEquivalentPlasticStrain(
      CDPMaterialTable::Branch::TENSION,
      linearized.result.state.tensile_equivalent_plastic_strain);
  const auto compression = _table.effectiveStrengthByEquivalentPlasticStrain(
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
  const double tension_step = 1.0e-7 * std::max(1.0, std::abs(tension.value));
  const double compression_step = 1.0e-7 * std::max(1.0, std::abs(compression.value));
  const double yield_tension_derivative =
      (yield_with_strengths(compression.value, tension.value + tension_step) -
       yield_with_strengths(compression.value, tension.value - tension_step)) /
      (2.0 * tension_step);
  const double yield_compression_derivative =
      (yield_with_strengths(compression.value + compression_step, tension.value) -
       yield_with_strengths(compression.value - compression_step, tension.value)) /
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
          yield_tension_derivative * tension.right_derivative / stress_scale;
    else
      direct[6] =
          yield_compression_derivative * compression.right_derivative / stress_scale;

    for (double & value : direct)
      value = -value;
    const auto unknown_derivative = solveLinearSystem(factorization, direct);
    ++linearized.result.local_backsolves;

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

AbaqusCDPLocalIntegrator::LocalJacobianDiagnostic
AbaqusCDPLocalIntegrator::localJacobianDiagnostic(const SymmetricTensor & total_strain,
                                                  const State & old_state) const
{
  const auto result = integrate(total_strain, old_state);
  if (!result.plastic)
    integrationError("local Jacobian diagnostic requires a plastic state");

  const double stress_scale = stressScale(total_strain, old_state);
  LocalVector unknown = {};
  for (std::size_t i = 0; i < 6; ++i)
    unknown[i] = result.effective_stress[i] / stress_scale;
  unknown[6] = result.plastic_multiplier / _strain_scale;
  unknown[7] = (result.state.tensile_equivalent_plastic_strain -
                old_state.tensile_equivalent_plastic_strain) /
               _strain_scale;
  unknown[8] = (result.state.compressive_equivalent_plastic_strain -
                old_state.compressive_equivalent_plastic_strain) /
               _strain_scale;
  return {automaticDifferentiationJacobian(unknown, total_strain, old_state, stress_scale),
          numericalJacobian(unknown, total_strain, old_state, stress_scale)};
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
