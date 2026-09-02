#pragma once

#include "AbaqusCDPLocalIntegrator.h"

/**
 * Transactional damage/recovery and Duvaut-Lions state layer for B-007/B-008.
 *
 * The inviscid backbone is advanced by AbaqusCDPLocalIntegrator. Plastic strain
 * and the two branch damage variables are then regularized independently. The
 * final Cauchy stress is assembled from the viscous effective stress and the
 * public w_t/w_c stiffness-recovery formula.
 */
class AbaqusCDPStateIntegrator
{
public:
  using SymmetricTensor = AbaqusCDPLocalIntegrator::SymmetricTensor;
  static constexpr std::size_t state_size = 16;
  static constexpr std::size_t transition_size = 6 + state_size;
  using TransitionColumn = std::array<double, transition_size>;
  using TransitionJacobian = std::array<TransitionColumn, transition_size>;

  struct Parameters
  {
    double tension_recovery;
    double compression_recovery;
    double relaxation_time;
    double state_tolerance = 1.0e-12;
  };

  struct State
  {
    AbaqusCDPLocalIntegrator::State backbone;
    SymmetricTensor viscous_plastic_strain = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double viscous_tension_damage = 0.0;
    double viscous_compression_damage = 0.0;
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
  };

  /** Rows are {Cauchy stress[6], new state[16]}; columns are
   * {total strain[6], old state[16]}. */
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

  /** Undamaged elastic stress associated with total_strain and backbone plastic history. */
  SymmetricTensor backboneEffectiveStress(const SymmetricTensor & total_strain,
                                          const State & state) const;

private:
  Result assembleResult(const SymmetricTensor & total_strain,
                        double time_step,
                        const State & old_state,
                        const AbaqusCDPLocalIntegrator::Result & backbone) const;

  const AbaqusCDPLocalIntegrator & _backbone_integrator;
  const Parameters _parameters;
};
