#include "CDPVelocityRateMaterial.h"
#include "ElasticVelocityRate.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include <algorithm>
#include <cmath>

registerMooseObject("DamSafetyApp", CDPVelocityRateMaterial);

InputParameters
CDPVelocityRateMaterial::validParams()
{
  auto params = Material::validParams();
  params.addClassDescription("Constant-C0 Newmark damping: 3D total rate or explicit 2D surrogate P; no CDP history mutation.");
  params.addRequiredCoupledVar("displacements", "Cartesian displacements");
  params.addRequiredCoupledVar("velocities", "Accepted Newmark velocity Aux variables");
  params.addRequiredCoupledVar("accelerations", "Accepted Newmark acceleration Aux variables");
  params.addCoupledVar("out_of_plane_strain", "Original static CDP weak-plane-stress unknown; never a damping-rate unknown");
  params.addParam<bool>("plane_stress_damping_projection", false, "Surrogate P: condense only the undamaged damping operator");
  params.addRequiredParam<Real>("beta", "Newmark beta used by inertia");
  params.addRequiredParam<Real>("gamma", "Newmark gamma used by inertia");
  params.addParam<Real>("alpha", 0, "HHT alpha used by inertia");
  params.addParam<MaterialPropertyName>("zeta", 0.0, "Constant stiffness damping coefficient");
  params.addRequiredParam<Real>("expected_youngs_modulus", "Undamaged isotropic Young modulus");
  params.addRequiredParam<Real>("expected_poissons_ratio", "Undamaged isotropic Poisson ratio");
  params.addRequiredParam<Real>("expected_stiffness_damping", "Constant nonnegative zeta");
  params.addParam<bool>("require_cdp_material", true, "Require CDP damage property; false only for elastic diagnostic controls");
  params.addParam<std::string>("base_name", "Material property prefix");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

CDPVelocityRateMaterial::CDPVelocityRateMaterial(const InputParameters & parameters)
  : Material(parameters),
    _prefix(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _beta(getParam<Real>("beta")), _gamma(getParam<Real>("gamma")),
    _alpha(getParam<Real>("alpha")), _youngs(getParam<Real>("expected_youngs_modulus")),
    _poisson(getParam<Real>("expected_poissons_ratio")),
    _expected_zeta(getParam<Real>("expected_stiffness_damping")),
    _plane_stress_projection(getParam<bool>("plane_stress_damping_projection")),
    _ndisp(_plane_stress_projection ? 2 : 3),
    _elasticity(getMaterialPropertyByName<RankFourTensor>(_prefix + "elasticity_tensor")),
    _zeta(getMaterialProperty<Real>("zeta")),
    _damage(getParam<bool>("require_cdp_material") ?
        &getMaterialPropertyByName<Real>(_prefix + "cdp_combined_damage") : nullptr),
    _thickness(_plane_stress_projection && isCoupled("out_of_plane_strain") ? &coupledValue("out_of_plane_strain") : nullptr),
    _thickness_old(_plane_stress_projection && isCoupled("out_of_plane_strain") ? &coupledValueOld("out_of_plane_strain") : nullptr),
    _rate(declareProperty<RankTwoTensor>(_prefix + (_plane_stress_projection ?
        "cdp_inplane_total_velocity_rate" : "cdp_total_velocity_rate"))),
    _damping_stress(declareProperty<RankTwoTensor>(_prefix + "cdp_velocity_damping_stress")),
    _projected_rate(declareProperty<RankTwoTensor>(_prefix + "cdp_damping_projected_velocity_rate")),
    _rate_tangent(declareProperty<RankFourTensor>(_prefix + "cdp_velocity_rate_tangent")),
    _thickness_increment_rate(declareProperty<Real>(_prefix + "cdp_static_thickness_increment_over_dt"))
{
  if (_mesh.dimension() != _ndisp || coupledComponents("displacements") != _ndisp ||
      coupledComponents("velocities") != _ndisp || coupledComponents("accelerations") != _ndisp ||
      (_plane_stress_projection && !_thickness) ||
      getParam<bool>("use_displaced_mesh") || !(_beta > 0) || !(_gamma > 0) ||
      !(_youngs > 0) || !(_poisson > -1 && _poisson < 0.5) ||
      !(_alpha >= -1.0 / 3 && _alpha <= 0) ||
      !std::isfinite(_beta) || !std::isfinite(_gamma) || !std::isfinite(_youngs) ||
      !std::isfinite(_expected_zeta) || _expected_zeta < 0)
    mooseError("CDP velocity diagnostic requires small-strain Cartesian Newmark data and valid constant coefficients");
  if (_plane_stress_projection)
  {
    auto * z = getVar("out_of_plane_strain", 0);
    if (z->kind() != Moose::VAR_SOLVER ||
        !((z->feType().family == LAGRANGE && z->feType().order == FIRST) ||
          (z->feType().family == MONOMIAL && z->feType().order == CONSTANT)))
      mooseError("CDP plane-stress projection requires FIRST LAGRANGE or CONSTANT MONOMIAL nonlinear thickness");
  }
  for (unsigned int c = 0; c < _ndisp; ++c)
  {
    _grad_u.push_back(&coupledGradient("displacements", c));
    _grad_u_old.push_back(&coupledGradientOld("displacements", c));
    _grad_v_old.push_back(&coupledGradientOld("velocities", c));
    _grad_a_old.push_back(&coupledGradientOld("accelerations", c));
  }
}

void
CDPVelocityRateMaterial::initialSetup()
{
  Material::initialSetup();
  if (_plane_stress_projection)
  {
    auto * z = getVar("out_of_plane_strain", 0);
    if (z->kind() != Moose::VAR_SOLVER ||
        !((z->feType().family == LAGRANGE && z->feType().order == FIRST) ||
          (z->feType().family == MONOMIAL && z->feType().order == CONSTANT)) || z->blockIDs() != getVar("displacements", 0)->blockIDs())
      mooseError("CDP plane-stress projection requires FIRST LAGRANGE or CONSTANT MONOMIAL nonlinear thickness in the displacement blocks");
  }
  for (unsigned int c = 0; c < _ndisp; ++c)
  {
    auto * u = getVar("displacements", c);
    auto * v = getVar("velocities", c);
    auto * a = getVar("accelerations", c);
    if (_plane_stress_projection && u->number() == getVar("out_of_plane_strain", 0)->number())
      mooseError("CDP plane-stress projection keeps static thickness distinct from in-plane displacements");
    if (u->kind() != Moose::VAR_SOLVER || v->kind() != Moose::VAR_AUXILIARY ||
        a->kind() != Moose::VAR_AUXILIARY || u->feType() != v->feType() ||
        u->feType() != a->feType() || u->blockIDs() != v->blockIDs() ||
        u->blockIDs() != a->blockIDs())
      mooseError("CDP velocity diagnostic requires accepted Aux histories in the displacement FE/block space");
  }
}

void
CDPVelocityRateMaterial::computeQpProperties()
{
  if (_coord_sys != Moose::COORD_XYZ)
    mooseError("CDP velocity damping requires Cartesian coordinates");
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
              std::abs(_elasticity[_qp](i,j,k,l) - expected) > 1e-10 * _youngs)
            mooseError("CDP velocity diagnostic requires the declared constant undamaged elasticity");
        }
  if (!std::isfinite(_zeta[_qp]) || _zeta[_qp] < 0 ||
      std::abs(_zeta[_qp] - _expected_zeta) > 1e-12 * std::max(1.0, _expected_zeta))
    mooseError("CDP velocity diagnostic requires the declared constant nonnegative damping coefficient");
  if (_damage && (!std::isfinite((*_damage)[_qp]) || (*_damage)[_qp] < 0 || (*_damage)[_qp] > 1))
    mooseError("CDP velocity diagnostic encountered invalid constitutive damage");
  RankTwoTensor gradient, old_gradient;
  for (unsigned int c = 0; c < _ndisp; ++c)
    for (unsigned int d = 0; d < _ndisp; ++d)
    {
      old_gradient(c,d) = (*_grad_v_old[c])[_qp](d);
      gradient(c,d) = _dt > 0 ? ElasticVelocityRate::trialVelocity(
          (*_grad_u[c])[_qp](d), (*_grad_u_old[c])[_qp](d),
          old_gradient(c,d), (*_grad_a_old[c])[_qp](d), _dt, _beta, _gamma)
          : old_gradient(c,d);
    }
  const RankTwoTensor old_rate = (old_gradient + old_gradient.transpose()) * 0.5;
  if (_t_step <= 1 && old_rate.L2norm() > 1e-14)
    mooseError("CDP velocity diagnostic requires zero initial strain rate");
  _rate[_qp] = (gradient + gradient.transpose()) * 0.5;
  _projected_rate[_qp] = _plane_stress_projection ?
      ElasticVelocityRate::closePlaneStressRate(_elasticity[_qp], _rate[_qp]) : _rate[_qp];
  const RankTwoTensor projected_old_rate = _plane_stress_projection ?
      ElasticVelocityRate::closePlaneStressRate(_elasticity[_qp], old_rate) : old_rate;
  _rate_tangent[_qp] = _plane_stress_projection ?
      ElasticVelocityRate::planeStressRateTangent(_elasticity[_qp]) : _elasticity[_qp];
  // This observable is a backward strain increment, NOT Newmark total thickness velocity.
  _thickness_increment_rate[_qp] = _thickness && _dt > 0 ?
      ((*_thickness)[_qp] - (*_thickness_old)[_qp]) / _dt : 0;
  // C0 stays undamaged. Jacobian_mult and damage/plastic increments are not damping inputs.
  _damping_stress[_qp] = _dt > 0 ? _elasticity[_qp] *
      (_projected_rate[_qp] * (1 + _alpha) - projected_old_rate * _alpha) * _zeta[_qp] : RankTwoTensor();
}
