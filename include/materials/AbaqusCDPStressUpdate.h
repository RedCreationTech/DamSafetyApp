#pragma once

#include "AbaqusCDPSubstepIntegrator.h"
#include "StressUpdateBase.h"
#include "CDPDiagnostics.h"
#include <fstream>

/**
 * MOOSE StressUpdateBase adapter for the table-driven Abaqus-CDP-compatible
 * material-point integrator.
 *
 * P0 is restricted to small-strain isotropic elasticity. All state changes are
 * reconstructed from old material properties on every nonlinear evaluation,
 * so a rejected global step cannot leak constitutive history.
 */
class AbaqusCDPStressUpdate : public StressUpdateBase
{
public:
  static InputParameters validParams();

  AbaqusCDPStressUpdate(const InputParameters & parameters);
  ~AbaqusCDPStressUpdate() override;

  void updateState(RankTwoTensor & strain_increment,
                   RankTwoTensor & inelastic_strain_increment,
                   const RankTwoTensor & rotation_increment,
                   RankTwoTensor & stress_new,
                   const RankTwoTensor & stress_old,
                   const RankFourTensor & elasticity_tensor,
                   const RankTwoTensor & elastic_strain_old,
                   bool compute_full_tangent_operator,
                   RankFourTensor & tangent_operator) override;

  bool requiresIsotropicTensor() override { return true; }
  bool isIsotropic() override { return true; }
  Real computeTimeStepLimit() override;
  TangentCalculationMethod getTangentCalculationMethod() override;

protected:
  void initQpStatefulProperties() override;
  void propagateQpStatefulProperties() override;

private:
  using CoreState = AbaqusCDPStateIntegrator::State;
  using SymmetricTensor = AbaqusCDPSubstepIntegrator::SymmetricTensor;

  static SymmetricTensor toArray(const RankTwoTensor & tensor);
  static RankTwoTensor toRankTwo(const SymmetricTensor & tensor);
  static void assignTangent(const AbaqusCDPSubstepIntegrator::TangentMatrix & source,
                            RankFourTensor & destination);
  CoreState oldState() const;
  void storeState(const AbaqusCDPSubstepIntegrator::LinearizedResult & result,
                  Real integration_microseconds);
  void writeCostSummary();
  void auditTangent(const SymmetricTensor & old_strain, const SymmetricTensor & new_strain,
                    const CoreState & old_state,
                    const AbaqusCDPSubstepIntegrator::LinearizedResult & result);

  const bool _enable_performance_diagnostics;
  const bool _enable_path_diagnostics;
  const bool _use_damage_t_secant_tangent;
  const Real _maximum_tensile_history_increment;
  const Real _minimum_state_timestep_limit;
  Real _state_timestep_limit;
  CDPDiagnostics::Counters _diagnostic_cost{};
  std::ofstream _diagnostic_trace;
  std::string _diagnostic_cost_file;
  unsigned int _failure_samples=0, _trace_calls=0, _tangent_checks=0;
  std::uint64_t _diagnostic_calls=0;
  const CDPMaterialTable _table;
  const AbaqusCDPLocalIntegrator _local_integrator;
  const AbaqusCDPStateIntegrator _state_integrator;
  const AbaqusCDPSubstepIntegrator _substep_integrator;

  MaterialProperty<RankTwoTensor> & _backbone_plastic_strain;
  const MaterialProperty<RankTwoTensor> & _backbone_plastic_strain_old;
  MaterialProperty<Real> & _tensile_equivalent_plastic_strain;
  const MaterialProperty<Real> & _tensile_equivalent_plastic_strain_old;
  MaterialProperty<Real> & _compressive_equivalent_plastic_strain;
  const MaterialProperty<Real> & _compressive_equivalent_plastic_strain_old;
  MaterialProperty<RankTwoTensor> & _viscous_plastic_strain;
  const MaterialProperty<RankTwoTensor> & _viscous_plastic_strain_old;
  MaterialProperty<Real> & _damage_t;
  const MaterialProperty<Real> & _damage_t_old;
  MaterialProperty<Real> & _damage_c;
  const MaterialProperty<Real> & _damage_c_old;

  MaterialProperty<Real> & _combined_damage;
  MaterialProperty<Real> & _stiffness_factor;
  MaterialProperty<Real> & _damage_t_secant_tangent_active;
  MaterialProperty<Real> & _local_iterations;
  MaterialProperty<Real> & _jacobian_fallbacks;
  MaterialProperty<Real> & _accepted_substeps;
  MaterialProperty<Real> & _failed_material_calls;
  MaterialProperty<Real> & _attempted_partitions;
  MaterialProperty<Real> & _maximum_partition_depth;
  MaterialProperty<Real> & _automatic_jacobian_evaluations;
  MaterialProperty<Real> & _finite_difference_jacobian_evaluations;
  MaterialProperty<Real> & _local_factorizations;
  MaterialProperty<Real> & _local_backsolves;
  MaterialProperty<Real> & _integration_microseconds;
};
