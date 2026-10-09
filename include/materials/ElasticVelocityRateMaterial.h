#pragma once
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Diagnostic-only elastic rate stresses. The constitutive S output is left untouched.
class ElasticVelocityRateMaterial : public Material
{
public:
  static InputParameters validParams();
  ElasticVelocityRateMaterial(const InputParameters & parameters);
  void initialSetup() override;
protected:
  void computeQpProperties() override;
  const std::string _prefix;
  const Real _beta, _gamma, _alpha, _youngs, _poisson;
  const bool _pointwise_rate;
  const MaterialProperty<RankFourTensor> & _elasticity;
  const MaterialProperty<RankFourTensor> & _tangent;
  const MaterialProperty<Real> & _zeta;
  const VariableValue * _rate_zz, * _rate_zz_old;
  std::vector<const VariableGradient *> _grad_u, _grad_u_old, _grad_v_old, _grad_a_old;
  MaterialProperty<RankTwoTensor> & _rate_stress, & _damping_stress;
  MaterialProperty<RankFourTensor> & _rate_tangent;
};
