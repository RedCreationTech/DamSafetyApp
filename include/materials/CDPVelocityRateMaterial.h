#pragma once
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Three-dimensional diagnostic: C0 times total strain rate, separate from constitutive S.
class CDPVelocityRateMaterial : public Material
{
public:
  static InputParameters validParams();
  CDPVelocityRateMaterial(const InputParameters & parameters);
  void initialSetup() override;
protected:
  void computeQpProperties() override;
  const std::string _prefix;
  const Real _beta, _gamma, _alpha, _youngs, _poisson, _expected_zeta;
  const MaterialProperty<RankFourTensor> & _elasticity;
  const MaterialProperty<Real> & _zeta;
  const MaterialProperty<Real> * _damage;
  std::vector<const VariableGradient *> _grad_u, _grad_u_old, _grad_v_old, _grad_a_old;
  MaterialProperty<RankTwoTensor> & _rate, & _damping_stress;
};
