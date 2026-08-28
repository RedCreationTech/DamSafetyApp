#pragma once

#include "AbaqusCDPStateIntegrator.h"

#include <array>

/**
 * Transactional binary substepping and numerical reference tangent for B-009A.
 *
 * This class differentiates the complete discrete material-point map. The
 * numerical matrix is a correctness reference, not the B-009B production
 * tangent intended for full-model performance.
 */
class AbaqusCDPSubstepIntegrator
{
public:
  using SymmetricTensor = AbaqusCDPStateIntegrator::SymmetricTensor;
  using State = AbaqusCDPStateIntegrator::State;
  using TangentMatrix = std::array<SymmetricTensor, 6>; // column-major

  struct Parameters
  {
    unsigned int maximum_substeps = 256;
    double maximum_strain_increment = 0.0;
    double tangent_perturbation = 1.0e-8;
  };

  struct Result
  {
    AbaqusCDPStateIntegrator::Result final_result;
    unsigned int accepted_substeps;
    unsigned int cutback_count;
    unsigned int attempted_partitions;
    unsigned int total_local_iterations;
    unsigned int total_jacobian_fallbacks;
    bool proactively_partitioned;
  };

  struct ReferenceTangent
  {
    TangentMatrix value;
    std::array<unsigned int, 6> plus_substeps;
    std::array<unsigned int, 6> minus_substeps;
    double perturbation;
  };

  struct LinearizedResult
  {
    Result result;
    TangentMatrix algorithmic_tangent;
  };

  AbaqusCDPSubstepIntegrator(const AbaqusCDPStateIntegrator & state_integrator,
                             Parameters parameters);

  Result integrate(const SymmetricTensor & old_total_strain,
                   const SymmetricTensor & new_total_strain,
                   double time_step,
                   const State & old_state) const;

  ReferenceTangent referenceTangent(const SymmetricTensor & old_total_strain,
                                    const SymmetricTensor & new_total_strain,
                                    double time_step,
                                    const State & old_state) const;

  LinearizedResult integrateLinearized(const SymmetricTensor & old_total_strain,
                                       const SymmetricTensor & new_total_strain,
                                       double time_step,
                                       const State & old_state) const;

  SymmetricTensor directionalDerivative(const SymmetricTensor & old_total_strain,
                                        const SymmetricTensor & new_total_strain,
                                        double time_step,
                                        const State & old_state,
                                        const SymmetricTensor & direction,
                                        double perturbation) const;

  static SymmetricTensor applyTangent(const TangentMatrix & tangent,
                                      const SymmetricTensor & direction);

private:
  const AbaqusCDPStateIntegrator & _state_integrator;
  const Parameters _parameters;
};
