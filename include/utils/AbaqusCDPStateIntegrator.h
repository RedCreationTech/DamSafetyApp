#pragma once

#include "AbaqusCDPLocalIntegrator.h"

/**
 * Transactional damage/recovery and Duvaut-Lions state layer for B-007/B-008.
 *
 * The inviscid backbone is advanced by AbaqusCDPLocalIntegrator. Plastic strain
 * and the two branch damage variables are regularized for branch output. An
 * opt-in compatibility candidate additionally regularizes the combined scalar
 * degradation used by the final Cauchy stress, matching the published Abaqus
 * Duvaut-Lions equation.
 */
class AbaqusCDPStateIntegrator
{
public:
  using SymmetricTensor = AbaqusCDPLocalIntegrator::SymmetricTensor;
  static constexpr std::size_t state_size = 17;
  static constexpr std::size_t transition_size = 6 + state_size;
  using TransitionColumn = std::array<double, transition_size>;
  using TransitionJacobian = std::array<TransitionColumn, transition_size>;

  struct Parameters
  {
    double tension_recovery;
    double compression_recovery;
    double relaxation_time;
    double state_tolerance = 1.0e-12;
    bool use_scalar_damage_viscosity = false;
  };

  struct State
  {
    AbaqusCDPLocalIntegrator::State backbone;
    SymmetricTensor viscous_plastic_strain = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double viscous_tension_damage = 0.0;
    double viscous_compression_damage = 0.0;
    double viscous_combined_damage = 0.0;
  };

  struct Result
  {
    AbaqusCDPLocalIntegrator::Result backbone;
    SymmetricTensor viscous_effective_stress;
    SymmetricTensor cauchy_stress;
    State state;
    AbaqusCDPFormula::DamageCombination damage;
    bool rate_independent;
    double dt_over_relaxation_time;
    double plastic_strain_lag_norm;
    double tension_damage_lag;
    double compression_damage_lag;
    double backbone_combined_damage;
    double combined_damage_lag;
  };

  /** Rows are {Cauchy stress[6], new state[17]}; columns are
   * {total strain[6], old state[17]}. */
  struct LinearizedResult
  {
    Result result;
    TransitionJacobian derivative;
  };

  AbaqusCDPStateIntegrator(const AbaqusCDPLocalIntegrator & backbone_integrator,
                           Parameters parameters);

  Result integrate(const SymmetricTensor & total_strain,
                   double time_step,
                   const State & old_state) const;
  LinearizedResult integrateLinearized(const SymmetricTensor & total_strain,
                                       double time_step,
                                       const State & old_state) const;

private:
  Result assembleResult(const SymmetricTensor & total_strain,
                        double time_step,
                        const State & old_state,
                        const AbaqusCDPLocalIntegrator::Result & backbone) const;

  const AbaqusCDPLocalIntegrator & _backbone_integrator;
  const Parameters _parameters;
};
