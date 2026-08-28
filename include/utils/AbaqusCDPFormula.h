#pragma once

#include <array>

/**
 * Stateless public-formula layer for an Abaqus-CDP-compatible material.
 *
 * Sign convention: tensile normal stress is positive and pressure
 * p = -trace(sigma) / 3 is positive in compression. Symmetric tensors use
 * {xx, yy, zz, xy, yz, xz}; shear entries are physical tensor components.
 * This layer does not perform return mapping or own constitutive state.
 */
namespace AbaqusCDPFormula
{
using SymmetricTensor = std::array<double, 6>;

struct StressInvariants
{
  double pressure;
  double mises;
  std::array<double, 3> principal_stress; // ascending
  double maximum_principal_stress;
  double tension_weight;
};

struct YieldCoefficients
{
  double alpha;
  double beta;
  double gamma;
};

struct FlowPotentialResult
{
  double value;
  /** Gradient for an ordinary dot product with {dxx,dyy,dzz,dxy,dyz,dxz}. */
  SymmetricTensor gradient;
};

struct DamageCombination
{
  double tension_weight;
  double tensile_recovery_factor;
  double compressive_recovery_factor;
  double stiffness_factor;
  double damage;
};

StressInvariants stressInvariants(const SymmetricTensor & effective_stress);

YieldCoefficients yieldCoefficients(double biaxial_to_uniaxial_compression_ratio,
                                    double tensile_meridian_ratio,
                                    double compression_strength,
                                    double tension_strength);

double yieldFunction(const SymmetricTensor & effective_stress,
                     double compression_strength,
                     double tension_strength,
                     double biaxial_to_uniaxial_compression_ratio,
                     double tensile_meridian_ratio);

FlowPotentialResult flowPotential(const SymmetricTensor & effective_stress,
                                  double dilation_angle_degrees,
                                  double eccentricity,
                                  double initial_tension_strength);

DamageCombination combineDamage(double compression_damage,
                                double tension_damage,
                                double tension_recovery,
                                double compression_recovery,
                                double tension_weight);

double duvautLionsUpdate(double old_viscous_state,
                         double backbone_state,
                         double time_step,
                         double relaxation_time);

SymmetricTensor duvautLionsUpdate(const SymmetricTensor & old_viscous_state,
                                  const SymmetricTensor & backbone_state,
                                  double time_step,
                                  double relaxation_time);
}
