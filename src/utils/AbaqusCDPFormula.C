#include "AbaqusCDPFormula.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace
{
constexpr double pi = 3.141592653589793238462643383279502884;

[[noreturn]] void
formulaError(const std::string & message)
{
  throw std::invalid_argument("AbaqusCDPFormula: " + message);
}

void
requireFinite(const double value, const std::string & name)
{
  if (!std::isfinite(value))
    formulaError(name + " must be finite");
}

double
positivePart(const double value)
{
  return std::max(0.0, value);
}

std::array<double, 3>
principalStress(const AbaqusCDPFormula::SymmetricTensor & stress)
{
  const double xx = stress[0];
  const double yy = stress[1];
  const double zz = stress[2];
  const double xy = stress[3];
  const double yz = stress[4];
  const double xz = stress[5];
  const double off_diagonal_norm = xy * xy + yz * yz + xz * xz;

  std::array<double, 3> result;
  if (off_diagonal_norm == 0.0)
    result = {xx, yy, zz};
  else
  {
    const double center = (xx + yy + zz) / 3.0;
    const double centered_norm = (xx - center) * (xx - center) +
                                 (yy - center) * (yy - center) +
                                 (zz - center) * (zz - center) + 2.0 * off_diagonal_norm;
    const double scale = std::sqrt(centered_norm / 6.0);

    if (scale == 0.0)
      result = {center, center, center};
    else
    {
      const double bxx = (xx - center) / scale;
      const double byy = (yy - center) / scale;
      const double bzz = (zz - center) / scale;
      const double bxy = xy / scale;
      const double byz = yz / scale;
      const double bxz = xz / scale;
      const double determinant = bxx * byy * bzz + 2.0 * bxy * byz * bxz -
                                 bxx * byz * byz - byy * bxz * bxz - bzz * bxy * bxy;
      const double angle = std::acos(std::clamp(determinant / 2.0, -1.0, 1.0)) / 3.0;
      const double largest = center + 2.0 * scale * std::cos(angle);
      const double smallest = center + 2.0 * scale * std::cos(angle + 2.0 * pi / 3.0);
      result = {smallest, 3.0 * center - largest - smallest, largest};
    }
  }

  std::sort(result.begin(), result.end());
  return result;
}

void
validateStress(const AbaqusCDPFormula::SymmetricTensor & stress)
{
  for (const double component : stress)
    requireFinite(component, "stress component");
}
}

namespace AbaqusCDPFormula
{
PrincipalSpectrum
stablePrincipalSpectrum(const SymmetricTensor & stress)
{
  validateStress(stress);
  using Matrix = std::array<std::array<long double, 3>, 3>;
  long double scale=0;
  for (auto x:stress) scale=std::max(scale,std::abs(static_cast<long double>(x)));
  if (scale==0) scale=1;
  Matrix a={{{stress[0]/scale,stress[3]/scale,stress[5]/scale},
             {stress[3]/scale,stress[1]/scale,stress[4]/scale},
             {stress[5]/scale,stress[4]/scale,stress[2]/scale}}};
  Matrix q={{{1,0,0},{0,1,0},{0,0,1}}};
  for (unsigned iteration=0;iteration<80;++iteration)
  {
    unsigned p=0,r=1;
    for (unsigned i=0;i<3;++i)
      for (unsigned j=i+1;j<3;++j)
        if (std::abs(a[i][j])>std::abs(a[p][r])) {p=i;r=j;}
    if (std::abs(a[p][r])<1e-30L) break;
    const long double tau=(a[r][r]-a[p][p])/(2*a[p][r]);
    const long double t=std::copysign(1.L,tau)/(std::abs(tau)+std::hypot(1.L,tau));
    const long double c=1/std::sqrt(1+t*t),s=t*c;
    const long double off=a[p][r];
    a[p][p]-=t*off; a[r][r]+=t*off; a[p][r]=a[r][p]=0;
    for (unsigned k=0;k<3;++k)
    {
      if (k!=p && k!=r)
      {
        const long double ap=a[k][p],ar=a[k][r];
        a[k][p]=a[p][k]=c*ap-s*ar; a[k][r]=a[r][k]=s*ap+c*ar;
      }
      const long double qp=q[k][p],qr=q[k][r];
      q[k][p]=c*qp-s*qr; q[k][r]=s*qp+c*qr;
    }
  }
  std::array<unsigned,3> order={{0,1,2}};
  std::sort(order.begin(),order.end(),[&](unsigned i,unsigned j){return a[i][i]<a[j][j];});
  PrincipalSpectrum result={};
  for(unsigned k=0;k<3;++k)
  {
    auto n=order[k]; result.values[k]=static_cast<double>(a[n][n]*scale);
    const long double x=q[0][n],y=q[1][n],z=q[2][n];
    result.gradients[k]={{static_cast<double>(x*x),static_cast<double>(y*y),static_cast<double>(z*z),
                         static_cast<double>(2*x*y),static_cast<double>(2*y*z),static_cast<double>(2*x*z)}};
  }
  return result;
}

StressInvariants
stressInvariants(const SymmetricTensor & effective_stress,
                 const bool use_stable_principal_stress)
{
  validateStress(effective_stress);
  const double trace = effective_stress[0] + effective_stress[1] + effective_stress[2];
  const double mean = trace / 3.0;
  const double sxx = effective_stress[0] - mean;
  const double syy = effective_stress[1] - mean;
  const double szz = effective_stress[2] - mean;
  const double j2 = 0.5 * (sxx * sxx + syy * syy + szz * szz +
                           2.0 * (effective_stress[3] * effective_stress[3] +
                                  effective_stress[4] * effective_stress[4] +
                                  effective_stress[5] * effective_stress[5]));
  const auto principal = use_stable_principal_stress
      ? stablePrincipalSpectrum(effective_stress).values : principalStress(effective_stress);

  double positive_sum = 0.0;
  double absolute_sum = 0.0;
  for (const double value : principal)
  {
    positive_sum += positivePart(value);
    absolute_sum += std::abs(value);
  }

  return {-mean,
          std::sqrt(std::max(0.0, 3.0 * j2)),
          principal,
          principal.back(),
          absolute_sum == 0.0 ? 0.0 : positive_sum / absolute_sum};
}

YieldCoefficients
yieldCoefficients(const double biaxial_to_uniaxial_compression_ratio,
                  const double tensile_meridian_ratio,
                  const double compression_strength,
                  const double tension_strength)
{
  requireFinite(biaxial_to_uniaxial_compression_ratio, "fb0/fc0");
  requireFinite(tensile_meridian_ratio, "Kc");
  requireFinite(compression_strength, "compression strength");
  requireFinite(tension_strength, "tension strength");
  if (biaxial_to_uniaxial_compression_ratio <= 0.5)
    formulaError("fb0/fc0 must be greater than 0.5");
  if (tensile_meridian_ratio <= 0.5 || tensile_meridian_ratio > 1.0)
    formulaError("Kc must be in (0.5, 1]");
  if (compression_strength <= 0.0 || tension_strength <= 0.0)
    formulaError("compression and tension strengths must be positive");

  const double alpha = (biaxial_to_uniaxial_compression_ratio - 1.0) /
                       (2.0 * biaxial_to_uniaxial_compression_ratio - 1.0);
  const double beta = compression_strength / tension_strength * (1.0 - alpha) -
                      (1.0 + alpha);
  const double gamma = 3.0 * (1.0 - tensile_meridian_ratio) /
                       (2.0 * tensile_meridian_ratio - 1.0);
  return {alpha, beta, gamma};
}

double
yieldFunction(const SymmetricTensor & effective_stress,
              const double compression_strength,
              const double tension_strength,
              const double biaxial_to_uniaxial_compression_ratio,
              const double tensile_meridian_ratio,
              const bool use_stable_principal_stress)
{
  const auto invariants = stressInvariants(effective_stress, use_stable_principal_stress);
  const auto coefficients = yieldCoefficients(biaxial_to_uniaxial_compression_ratio,
                                               tensile_meridian_ratio,
                                               compression_strength,
                                               tension_strength);
  const double maximum = invariants.maximum_principal_stress;
  return (invariants.mises - 3.0 * coefficients.alpha * invariants.pressure +
          coefficients.beta * positivePart(maximum) -
          coefficients.gamma * positivePart(-maximum)) /
             (1.0 - coefficients.alpha) -
         compression_strength;
}

FlowPotentialResult
flowPotential(const SymmetricTensor & effective_stress,
              const double dilation_angle_degrees,
              const double eccentricity,
              const double initial_tension_strength)
{
  validateStress(effective_stress);
  requireFinite(dilation_angle_degrees, "dilation angle");
  requireFinite(eccentricity, "eccentricity");
  requireFinite(initial_tension_strength, "initial tension strength");
  if (dilation_angle_degrees < 0.0 || dilation_angle_degrees >= 90.0)
    formulaError("dilation angle must be in [0, 90) degrees");
  if (eccentricity <= 0.0 || initial_tension_strength <= 0.0)
    formulaError("eccentricity and initial tension strength must be positive");

  const auto invariants = stressInvariants(effective_stress);
  const double tangent = std::tan(dilation_angle_degrees * pi / 180.0);
  const double asymptote = eccentricity * initial_tension_strength * tangent;
  const double root = std::sqrt(asymptote * asymptote + invariants.mises * invariants.mises);

  SymmetricTensor gradient = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  if (root > 0.0)
  {
    const double mean = (effective_stress[0] + effective_stress[1] + effective_stress[2]) / 3.0;
    const double normal_factor = 1.5 / root;
    gradient[0] = normal_factor * (effective_stress[0] - mean);
    gradient[1] = normal_factor * (effective_stress[1] - mean);
    gradient[2] = normal_factor * (effective_stress[2] - mean);
    gradient[3] = 3.0 / root * effective_stress[3];
    gradient[4] = 3.0 / root * effective_stress[4];
    gradient[5] = 3.0 / root * effective_stress[5];
  }
  gradient[0] += tangent / 3.0;
  gradient[1] += tangent / 3.0;
  gradient[2] += tangent / 3.0;

  return {root - invariants.pressure * tangent, gradient};
}

DamageCombination
combineDamage(const double compression_damage,
              const double tension_damage,
              const double tension_recovery,
              const double compression_recovery,
              const double tension_weight)
{
  requireFinite(compression_damage, "compression damage");
  requireFinite(tension_damage, "tension damage");
  requireFinite(tension_recovery, "tension recovery");
  requireFinite(compression_recovery, "compression recovery");
  requireFinite(tension_weight, "tension weight");
  if (compression_damage < 0.0 || compression_damage > 1.0 || tension_damage < 0.0 ||
      tension_damage > 1.0)
    formulaError("damage values must be in [0, 1]");
  if (tension_recovery < 0.0 || tension_recovery > 1.0 || compression_recovery < 0.0 ||
      compression_recovery > 1.0)
    formulaError("recovery factors must be in [0, 1]");
  if (tension_weight < 0.0 || tension_weight > 1.0)
    formulaError("tension weight must be in [0, 1]");

  const double tensile_recovery_factor = 1.0 - tension_recovery * tension_weight;
  const double compressive_recovery_factor =
      1.0 - compression_recovery * (1.0 - tension_weight);
  const double stiffness_factor = (1.0 - tensile_recovery_factor * compression_damage) *
                                  (1.0 - compressive_recovery_factor * tension_damage);
  return {tension_weight,
          tensile_recovery_factor,
          compressive_recovery_factor,
          stiffness_factor,
          1.0 - stiffness_factor};
}

double
duvautLionsUpdate(const double old_viscous_state,
                  const double backbone_state,
                  const double time_step,
                  const double relaxation_time)
{
  requireFinite(old_viscous_state, "old viscous state");
  requireFinite(backbone_state, "backbone state");
  requireFinite(time_step, "time step");
  requireFinite(relaxation_time, "relaxation time");
  if (time_step < 0.0 || relaxation_time < 0.0)
    formulaError("time step and relaxation time must be nonnegative");
  if (relaxation_time == 0.0)
    return backbone_state;
  if (time_step == 0.0)
    return old_viscous_state;
  return (relaxation_time * old_viscous_state + time_step * backbone_state) /
         (relaxation_time + time_step);
}

SymmetricTensor
duvautLionsUpdate(const SymmetricTensor & old_viscous_state,
                  const SymmetricTensor & backbone_state,
                  const double time_step,
                  const double relaxation_time)
{
  SymmetricTensor result;
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = duvautLionsUpdate(
        old_viscous_state[i], backbone_state[i], time_step, relaxation_time);
  return result;
}
}
