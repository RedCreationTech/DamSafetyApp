#pragma once
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// C0 damping: 3D total rate or 2D surrogate P, separate from constitutive S.
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
  const bool _plane_stress_projection;
  const unsigned int _ndisp;
  const MaterialProperty<RankFourTensor> & _elasticity;
  const MaterialProperty<Real> & _zeta;
  const MaterialProperty<Real> * _damage;
  const VariableValue * _thickness, * _thickness_old;
  std::vector<const VariableGradient *> _grad_u, _grad_u_old, _grad_v_old, _grad_a_old;
  MaterialProperty<RankTwoTensor> & _rate, & _damping_stress;
  MaterialProperty<RankTwoTensor> & _projected_rate;
  MaterialProperty<RankFourTensor> & _rate_tangent;
  MaterialProperty<Real> & _thickness_increment_rate;
};
