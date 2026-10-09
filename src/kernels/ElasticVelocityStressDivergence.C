#include "ElasticVelocityStressDivergence.h"
#include "ElasticVelocityRate.h"
#include "ElasticityTensorTools.h"

registerMooseObject("DamSafetyApp", ElasticVelocityStressDivergence);

InputParameters
ElasticVelocityStressDivergence::validParams()
{
  auto params = DynamicStressDivergenceTensors::validParams();
  params.addClassDescription("Strict elastic HHT stress divergence using reconstructed Newmark velocity; constitutive stress output is unchanged.");
  params.addRequiredParam<Real>("beta", "Same Newmark beta as inertia");
  params.addRequiredParam<Real>("gamma", "Same Newmark gamma as inertia");
  params.addCoupledVar("out_of_plane_strain_rate", "Algebraic weak thickness rate; absent for audited observable elimination");
  return params;
}

ElasticVelocityStressDivergence::ElasticVelocityStressDivergence(const InputParameters & parameters)
  : DynamicStressDivergenceTensors(parameters),
    _beta(getParam<Real>("beta")), _gamma(getParam<Real>("gamma")),
    _rate_var(isCoupled("out_of_plane_strain_rate") ? coupled("out_of_plane_strain_rate") : 0),
    _rate_coupled(isCoupled("out_of_plane_strain_rate")),
    _elasticity(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
    _rate_tangent(getMaterialPropertyByName<RankFourTensor>(_base_name + "elastic_velocity_rate_tangent")),
    _damping_stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "elastic_velocity_damping_stress"))
{
  if (_ndisp != 2 || _component > 1 || !_out_of_plane_strain_coupled ||
      _out_of_plane_direction != 2 || _use_finite_deform_jacobian ||
      _volumetric_locking_correction || _use_displaced_mesh || _static_initialization ||
      isCoupled("temperature") || !getParam<std::vector<MaterialPropertyName>>("eigenstrain_names").empty() ||
      !(_beta > 0) || !(_gamma > 0))
    mooseError("Elastic velocity divergence supports only incremental SMALL XY weak plane stress without Bbar, eigenstrain, temperature or static initialization");
}

Real
ElasticVelocityStressDivergence::computeQpResidual()
{
  if (_dt <= 0)
    return 0;
  const auto stress = _stress[_qp] * (1 + _alpha) - _stress_old[_qp] * _alpha + _damping_stress[_qp];
  return stress.row(_component) * _grad_test[_i][_qp];
}

Real
ElasticVelocityStressDivergence::computeQpJacobian()
{
  if (_dt <= 0)
    return 0;
  return (1 + _alpha) * (StressDivergenceTensors::computeQpJacobian() +
      _zeta[_qp] * ElasticVelocityRate::displacementDerivative(_dt, _beta, _gamma) *
      ElasticityTensorTools::elasticJacobian(_rate_tangent[_qp], _component, _component,
                                            _grad_test[_i][_qp], _grad_phi[_j][_qp]));
}

Real
ElasticVelocityStressDivergence::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_dt <= 0)
    return 0;
  if (_rate_coupled && jvar == _rate_var)
  {
    Real derivative = 0;
    for (unsigned int d = 0; d < 2; ++d)
      derivative += _elasticity[_qp](_component,d,2,2) * _grad_test[_i][_qp](d);
    return _zeta[_qp] * (1 + _alpha) * derivative * _phi[_j][_qp];
  }
  for (unsigned int c = 0; c < _ndisp; ++c)
    if (jvar == _disp_var[c])
      return (1 + _alpha) * (StressDivergenceTensors::computeQpOffDiagJacobian(jvar) +
          _zeta[_qp] * ElasticVelocityRate::displacementDerivative(_dt, _beta, _gamma) *
          ElasticityTensorTools::elasticJacobian(_rate_tangent[_qp], _component, c,
                                                _grad_test[_i][_qp], _grad_phi[_j][_qp]));
  // Original thickness strain belongs only to static stress, not trial velocity.
  return (1 + _alpha) * StressDivergenceTensors::computeQpOffDiagJacobian(jvar);
}
