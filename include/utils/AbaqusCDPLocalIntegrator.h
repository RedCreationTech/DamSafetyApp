#pragma once

#include "AbaqusCDPFormula.h"
#include "CDPMaterialTable.h"

#include <string>

/**
 * Material-point effective-stress integrator for the B-006B prototype.
 *
 * The class is intentionally independent of MOOSE material state. A caller
 * passes an immutable old state and receives a new state only after local
 * convergence, which makes cutback/rollback transactional by construction.
 * Damage, viscosity history, and the global algorithmic tangent are outside
 * this B-006B layer.
 */
class AbaqusCDPLocalIntegrator
{
public:
  using SymmetricTensor = AbaqusCDPFormula::SymmetricTensor;

  enum class ActiveBranch
  {
    ELASTIC,
    TENSION,
    COMPRESSION,
    MIXED
  };

  struct Parameters
  {
    double youngs_modulus;
    double poissons_ratio;
    double dilation_angle_degrees;
    double eccentricity;
    double biaxial_to_uniaxial_compression_ratio;
    double tensile_meridian_ratio;
    unsigned int maximum_iterations = 40;
    double residual_tolerance = 1.0e-9;
    double finite_difference_step = 1.0e-7;
    double minimum_line_search = 1.0e-6;
  };

  struct State
  {
    SymmetricTensor plastic_strain = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double tensile_equivalent_plastic_strain = 0.0;
    double compressive_equivalent_plastic_strain = 0.0;
  };

  struct Result
  {
    SymmetricTensor effective_stress;
    State state;
    ActiveBranch active_branch;
    bool plastic;
    unsigned int iterations;
    double residual_norm;
    double trial_yield;
    double final_yield;
    double plastic_multiplier;
    double backbone_tension_damage;
    double backbone_compression_damage;
  };

  AbaqusCDPLocalIntegrator(const CDPMaterialTable & table, Parameters parameters);

  Result integrate(const SymmetricTensor & total_strain, const State & old_state) const;
  SymmetricTensor elasticStress(const SymmetricTensor & elastic_strain) const;

  static std::string branchName(ActiveBranch branch);

private:
  const CDPMaterialTable & _table;
  const Parameters _parameters;
  const double _initial_tension_strength;
  const double _initial_compression_strength;
  const double _stress_scale_floor;
  const double _strain_scale;
};
