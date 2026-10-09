#include "ElasticPlaneStressRate.h"
#include "ElasticVelocityRate.h"

registerMooseObject("DamSafetyApp", ElasticPlaneStressRate);

InputParameters
ElasticPlaneStressRate::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription("Algebraic thickness strain rate constraint in the original weak-plane-stress space; no thickness inertia.");
  params.addRequiredCoupledVar("displacements", "In-plane displacements");
  params.addRequiredParam<Real>("beta", "Same Newmark beta as inertia");
  params.addRequiredParam<Real>("gamma", "Same Newmark gamma as inertia");
  params.addParam<std::string>("base_name", "Material property prefix");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

ElasticPlaneStressRate::ElasticPlaneStressRate(const InputParameters & parameters)
  : Kernel(parameters), _beta(getParam<Real>("beta")), _gamma(getParam<Real>("gamma")),
    _elasticity(getMaterialPropertyByName<RankFourTensor>(
        (isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "") + "elasticity_tensor")),
    _rate_stress(getMaterialPropertyByName<RankTwoTensor>(
        (isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "") + "elastic_velocity_rate_stress"))
{
  if (coupledComponents("displacements") != 2 || getParam<bool>("use_displaced_mesh") ||
      !(_beta > 0) || !(_gamma > 0))
    mooseError("Elastic plane stress rate requires two in-plane components and valid Newmark constants");
  for (unsigned int c = 0; c < 2; ++c)
    _displacements.push_back(coupled("displacements", c));
}

Real
ElasticPlaneStressRate::computeQpResidual()
{
  return _rate_stress[_qp](2,2) * _test[_i][_qp];
}

Real
ElasticPlaneStressRate::computeQpJacobian()
{
  return _elasticity[_qp](2,2,2,2) * _phi[_j][_qp] * _test[_i][_qp];
}

Real
ElasticPlaneStressRate::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_dt <= 0)
    return 0;
  for (unsigned int c = 0; c < 2; ++c)
    if (jvar == _displacements[c])
    {
      Real derivative = 0;
      for (unsigned int d = 0; d < 2; ++d)
        derivative += _elasticity[_qp](2,2,c,d) * _grad_phi[_j][_qp](d);
      return _test[_i][_qp] * derivative *
          ElasticVelocityRate::displacementDerivative(_dt, _beta, _gamma);
    }
  return 0;
}
