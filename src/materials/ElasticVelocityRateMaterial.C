#include "ElasticVelocityRateMaterial.h"
#include "ElasticVelocityRate.h"
#include "MooseVariable.h"
#include "MooseMesh.h"
#include <algorithm>
#include <cmath>

registerMooseObject("DamSafetyApp", ElasticVelocityRateMaterial);

InputParameters
ElasticVelocityRateMaterial::validParams()
{
  auto params = Material::validParams();
  params.addClassDescription("Strictly elastic Newmark trial velocity rate and HHT damping stress; no constitutive history is changed.");
  params.addRequiredCoupledVar("displacements", "Two in-plane displacements");
  params.addRequiredCoupledVar("velocities", "Accepted Newmark velocity Aux variables");
  params.addRequiredCoupledVar("accelerations", "Accepted Newmark acceleration Aux variables");
  params.addRequiredCoupledVar("out_of_plane_strain", "Original weak-plane-stress unknown");
  params.addRequiredCoupledVar("out_of_plane_strain_rate", "Algebraic thickness rate in the same FE/block space");
  params.addRequiredParam<Real>("beta", "Newmark beta, identical to inertia/Aux update");
  params.addRequiredParam<Real>("gamma", "Newmark gamma, identical to inertia/Aux update");
  params.addParam<Real>("alpha", 0, "HHT alpha, identical to inertia");
  params.addParam<MaterialPropertyName>("zeta", 0.0, "Elastic stiffness-proportional damping coefficient");
  params.addRequiredParam<Real>("expected_youngs_modulus", "Audited constant elastic Young modulus");
  params.addRequiredParam<Real>("expected_poissons_ratio", "Audited constant elastic Poisson ratio");
  params.addRequiredParam<bool>("strict_elastic_model", "Explicit acknowledgement: constant isotropic elasticity, incremental small strain, no damage/eigenstrain");
  params.addParam<std::string>("base_name", "Material property prefix");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

ElasticVelocityRateMaterial::ElasticVelocityRateMaterial(const InputParameters & parameters)
  : Material(parameters),
    _prefix(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _beta(getParam<Real>("beta")), _gamma(getParam<Real>("gamma")),
    _alpha(getParam<Real>("alpha")), _youngs(getParam<Real>("expected_youngs_modulus")),
    _poisson(getParam<Real>("expected_poissons_ratio")),
    _elasticity(getMaterialPropertyByName<RankFourTensor>(_prefix + "elasticity_tensor")),
    _tangent(getMaterialPropertyByName<RankFourTensor>(_prefix + "Jacobian_mult")),
    _zeta(getMaterialProperty<Real>("zeta")),
    _rate_zz(coupledValue("out_of_plane_strain_rate")),
    _rate_zz_old(coupledValueOld("out_of_plane_strain_rate")),
    _rate_stress(declareProperty<RankTwoTensor>(_prefix + "elastic_velocity_rate_stress")),
    _damping_stress(declareProperty<RankTwoTensor>(_prefix + "elastic_velocity_damping_stress"))
{
  if (!getParam<bool>("strict_elastic_model") || getParam<bool>("use_displaced_mesh") ||
      !(_beta > 0) || !(_gamma > 0) || !(_youngs > 0) ||
      !(_poisson > -1 && _poisson < 0.5) || !(_alpha >= -1.0 / 3 && _alpha <= 0))
    mooseError("ElasticVelocityRateMaterial requires acknowledged constant small-strain elasticity and valid integration constants");
  if (coupledComponents("displacements") != 2 || coupledComponents("velocities") != 2 ||
      coupledComponents("accelerations") != 2 || _mesh.dimension() != 2)
    mooseError("ElasticVelocityRateMaterial requires exactly two displacement/velocity/acceleration components");
  for (unsigned int c = 0; c < 2; ++c)
  {
    _grad_u.push_back(&coupledGradient("displacements", c));
    _grad_u_old.push_back(&coupledGradientOld("displacements", c));
    _grad_v_old.push_back(&coupledGradientOld("velocities", c));
    _grad_a_old.push_back(&coupledGradientOld("accelerations", c));
  }
}

void
ElasticVelocityRateMaterial::initialSetup()
{
  Material::initialSetup();
  auto * z = getVar("out_of_plane_strain", 0);
  auto * r = getVar("out_of_plane_strain_rate", 0);
  if (z->number() == r->number() || z->feType() != r->feType() || z->blockIDs() != r->blockIDs())
    mooseError("Elastic thickness strain and rate must be distinct nonlinear variables in the same FE/block space");
  if (z->kind() != Moose::VAR_SOLVER || r->kind() != Moose::VAR_SOLVER)
    mooseError("Elastic thickness strain and rate must be nonlinear variables");
  for (unsigned int c = 0; c < 2; ++c)
  {
    auto * u = getVar("displacements", c);
    if (u->number() == r->number() || u->number() == z->number() ||
        getVar("velocities", c)->kind() != Moose::VAR_AUXILIARY ||
        getVar("accelerations", c)->kind() != Moose::VAR_AUXILIARY ||
        u->feType() != getVar("velocities", c)->feType() ||
        u->feType() != getVar("accelerations", c)->feType() ||
        u->blockIDs() != getVar("velocities", c)->blockIDs() ||
        u->blockIDs() != getVar("accelerations", c)->blockIDs())
      mooseError("Newmark history must use the displacement FE/block space and distinct thickness unknowns");
  }
  if (hasMaterialPropertyByName<Real>(_prefix + "DamageT") ||
      hasMaterialPropertyByName<Real>(_prefix + "DamageC") ||
      hasMaterialPropertyByName<Real>(_prefix + "cdp_combined_damage"))
    mooseError("Elastic velocity damping candidate does not support CDP/damage");
}

void
ElasticVelocityRateMaterial::computeQpProperties()
{
  // Checking both tensors prevents accidentally using a degraded/rotated tangent or a
  // nonconstant elastic tensor. This candidate is deliberately restricted, not a CDP model.
  const Real mu = _youngs / (2 * (1 + _poisson));
  const Real lambda = _youngs * _poisson / ((1 + _poisson) * (1 - 2 * _poisson));
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int l = 0; l < 3; ++l)
        {
          const Real expected = lambda * (i == j) * (k == l) +
              mu * ((i == k && j == l) + (i == l && j == k));
          if (!std::isfinite(_elasticity[_qp](i,j,k,l)) ||
              !std::isfinite(_tangent[_qp](i,j,k,l)) ||
              std::abs(_elasticity[_qp](i,j,k,l) - expected) > 1e-10 * _youngs ||
              std::abs(_tangent[_qp](i,j,k,l) - expected) > 1e-10 * _youngs)
            mooseError("Elastic velocity damping requires the declared constant isotropic elasticity and tangent");
        }
  if (!std::isfinite(_zeta[_qp]) || _zeta[_qp] < 0)
    mooseError("Elastic velocity damping requires finite nonnegative zeta");
  RankTwoTensor gradient, old_gradient;
  for (unsigned int c = 0; c < 2; ++c)
    for (unsigned int d = 0; d < 2; ++d)
    {
      old_gradient(c,d) = (*_grad_v_old[c])[_qp](d);
      gradient(c,d) = _dt > 0 ? ElasticVelocityRate::trialVelocity(
          (*_grad_u[c])[_qp](d), (*_grad_u_old[c])[_qp](d),
          old_gradient(c,d), (*_grad_a_old[c])[_qp](d), _dt, _beta, _gamma)
          : old_gradient(c,d);
    }
  const RankTwoTensor old_rate = (old_gradient + old_gradient.transpose()) * 0.5;
  if (_t_step <= 1 && (old_rate.L2norm() > 1e-14 || std::abs(_rate_zz_old[_qp]) > 1e-14))
    mooseError("Elastic velocity candidate requires zero initial strain rate; nonzero initial acceleration is retained");
  RankTwoTensor rate = (gradient + gradient.transpose()) * 0.5;
  rate(2,2) = _rate_zz[_qp];
  RankTwoTensor previous_rate = old_rate;
  previous_rate(2,2) = _rate_zz_old[_qp];
  _rate_stress[_qp] = _elasticity[_qp] * rate;
  _damping_stress[_qp] = _dt > 0 ? _elasticity[_qp] *
      (rate * (1 + _alpha) - previous_rate * _alpha) * _zeta[_qp] : RankTwoTensor();
}
