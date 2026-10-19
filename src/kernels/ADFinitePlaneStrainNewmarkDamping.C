#include "ADFinitePlaneStrainNewmarkDamping.h"
#include "FinitePlaneStrainVelocityRate.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "MooseException.h"
#include <algorithm>
#include <cmath>

registerMooseObject("DamSafetyApp", ADFinitePlaneStrainNewmarkDamping);

InputParameters
ADFinitePlaneStrainNewmarkDamping::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Opt-in finite plane-strain constant-C0 Newmark damping, reference weak form; no CDP state mutation.");
  params.addRequiredParam<unsigned int>("component", "In-plane displacement component");
  params.addRequiredCoupledVar("displacements", "Two in-plane displacements");
  params.addRequiredCoupledVar("velocities", "Accepted Newmark auxiliary velocities");
  params.addRequiredCoupledVar("accelerations", "Accepted Newmark auxiliary accelerations");
  params.addRequiredParam<Real>("beta", "Newmark beta");
  params.addRequiredParam<Real>("gamma", "Newmark gamma");
  params.addRequiredParam<Real>("alpha", "HHT alpha");
  params.addRequiredParam<Real>("expected_youngs_modulus", "Constant C0 Young modulus");
  params.addRequiredParam<Real>("expected_poissons_ratio", "Constant C0 Poisson ratio");
  params.addRequiredParam<Real>("expected_stiffness_damping", "Constant nonnegative original stiffness damping");
  params.addRequiredParam<MaterialPropertyName>("zeta", "Original stiffness damping coefficient or property");
  params.addParam<MaterialPropertyName>("elasticity_tensor", "elasticity_tensor", "Undamaged isotropic C0");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

ADFinitePlaneStrainNewmarkDamping::ADFinitePlaneStrainNewmarkDamping(const InputParameters & parameters)
  : ADKernel(parameters), _component(getParam<unsigned int>("component")),
    _beta(getParam<Real>("beta")), _gamma(getParam<Real>("gamma")),
    _alpha(getParam<Real>("alpha")), _youngs(getParam<Real>("expected_youngs_modulus")),
    _poisson(getParam<Real>("expected_poissons_ratio")),
    _expected_zeta(getParam<Real>("expected_stiffness_damping")),
    _zeta(getMaterialProperty<Real>("zeta")),
    _elasticity(getMaterialProperty<RankFourTensor>("elasticity_tensor"))
{
  if (_mesh.dimension() != 2 || _component >= 2 || _use_displaced_mesh ||
      coupledComponents("displacements") != 2 || coupledComponents("velocities") != 2 ||
      coupledComponents("accelerations") != 2 ||
      !(_beta > 0) || !(_gamma > 0) || !(_alpha >= -1.0/3 && _alpha <= 0) ||
      !(_youngs > 0) || !(_poisson > -1 && _poisson < .5) ||
      !std::isfinite(_beta) || !std::isfinite(_gamma) || !std::isfinite(_youngs) ||
      !std::isfinite(_expected_zeta) || _expected_zeta < 0)
    mooseError("Finite Newmark damping requires 2D XY plane strain, reference mesh and finite valid constant coefficients");
  if (_var.number() != coupled("displacements", _component))
    paramError("variable", "Must match the selected displacement component");
  for (unsigned int c = 0; c < 2; ++c)
  {
    _gradients.push_back(&adCoupledGradient("displacements", c));
    _old_gradients.push_back(&coupledGradientOld("displacements", c));
    _old_velocities.push_back(&coupledGradientOld("velocities", c));
    _old_accelerations.push_back(&coupledGradientOld("accelerations", c));
  }
}

void
ADFinitePlaneStrainNewmarkDamping::initialSetup()
{
  if (getBlockCoordSystem() != Moose::COORD_XYZ)
    mooseError("Finite Newmark damping requires Cartesian coordinates");
  for (unsigned int c = 0; c < 2; ++c)
  {
    auto * u = getVar("displacements", c);
    auto * v = getVar("velocities", c);
    auto * a = getVar("accelerations", c);
    if (u->kind() != Moose::VAR_SOLVER || v->kind() != Moose::VAR_AUXILIARY ||
        a->kind() != Moose::VAR_AUXILIARY || u->feType() != v->feType() ||
        u->feType() != a->feType() || u->blockIDs() != v->blockIDs() ||
        u->blockIDs() != a->blockIDs())
      mooseError("Finite Newmark damping requires matching displacement and auxiliary FE/block spaces");
  }
}

ADReal
ADFinitePlaneStrainNewmarkDamping::computeQpResidual()
{
  if (!std::isfinite(_zeta[_qp]) ||
      std::abs(_zeta[_qp] - _expected_zeta) > 1e-12 * std::max(1.0, _expected_zeta))
    mooseError("Finite Newmark damping requires the declared constant stiffness damping");
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
            mooseError("Finite Newmark damping requires the declared constant undamaged isotropic C0");
        }
  ADRankTwoTensor gradient;
  RankTwoTensor old_gradient, old_velocity, old_acceleration;
  for (unsigned int c = 0; c < 2; ++c)
    for (unsigned int d = 0; d < 2; ++d)
    {
      gradient(c,d) = (*_gradients[c])[_qp](d);
      old_gradient(c,d) = (*_old_gradients[c])[_qp](d);
      old_velocity(c,d) = (*_old_velocities[c])[_qp](d);
      old_acceleration(c,d) = (*_old_accelerations[c])[_qp](d);
    }
  try
  {
    const auto nominal = FinitePlaneStrainVelocityRate::nominalStress(
        gradient, old_gradient, old_velocity, old_acceleration, _dt,
        _beta, _gamma, _alpha, _zeta[_qp], _youngs, _poisson);
    return nominal.row(_component) * _grad_test[_i][_qp];
  }
  catch (const std::invalid_argument & error)
  {
    throw MooseException(error.what());
  }
}
