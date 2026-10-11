#include "CDPVelocityStressDivergence.h"
#include "ElasticVelocityRate.h"
#include "ElasticityTensorTools.h"
#include "MooseMesh.h"

registerMooseObject("DamSafetyApp", CDPVelocityStressDivergence);

InputParameters
CDPVelocityStressDivergence::validParams()
{
  auto params = DynamicStressDivergenceTensors::validParams();
  params.addClassDescription("CDP HHT balance with 3D C0 total rate or 2D projected damping; static thickness balance is unchanged.");
  params.addRequiredParam<Real>("beta", "Newmark beta used by inertia");
  params.addRequiredParam<Real>("gamma", "Newmark gamma used by inertia");
  params.addParam<bool>("plane_stress_damping_projection", false, "2D surrogate P, without changing static CDP thickness balance");
  return params;
}

CDPVelocityStressDivergence::CDPVelocityStressDivergence(const InputParameters & parameters)
  : DynamicStressDivergenceTensors(parameters),
    _beta(getParam<Real>("beta")), _gamma(getParam<Real>("gamma")),
    _plane_stress_projection(getParam<bool>("plane_stress_damping_projection")),
    _elasticity(getMaterialPropertyByName<RankFourTensor>(_base_name + "cdp_velocity_rate_tangent")),
    _damping_stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "cdp_velocity_damping_stress"))
{
  const bool supported_dimension = _plane_stress_projection ?
      (_ndisp == 2 && _mesh.dimension() == 2 && _component < 2 &&
       _out_of_plane_strain_coupled && _out_of_plane_direction == 2) :
      (_ndisp == 3 && _mesh.dimension() == 3 && _component < 3 && !_out_of_plane_strain_coupled);
  if (!supported_dimension || _use_finite_deform_jacobian ||
      _volumetric_locking_correction || _use_displaced_mesh || _static_initialization ||
      isCoupled("temperature") ||
      !getParam<std::vector<MaterialPropertyName>>("eigenstrain_names").empty() ||
      !(_beta > 0) || !(_gamma > 0))
    mooseError("CDP velocity divergence supports incremental SMALL 3D or XY weak-plane-stress projection without Bbar, eigenstrain, temperature or static initialization");
}

Real
CDPVelocityStressDivergence::computeQpResidual()
{
  // Preserve the original arithmetic and nonlinear path at exact zero damping.
  if (_zeta[_qp] == 0.0)
    return DynamicStressDivergenceTensors::computeQpResidual();
  if (_dt <= 0)
    return 0;
  const auto stress = _stress[_qp] * (1 + _alpha) - _stress_old[_qp] * _alpha + _damping_stress[_qp];
  return stress.row(_component) * _grad_test[_i][_qp];
}

Real
CDPVelocityStressDivergence::computeQpJacobian()
{
  if (_zeta[_qp] == 0.0)
    return DynamicStressDivergenceTensors::computeQpJacobian();
  if (_dt <= 0)
    return 0;
  return (1 + _alpha) * (StressDivergenceTensors::computeQpJacobian() +
      _zeta[_qp] * ElasticVelocityRate::displacementDerivative(_dt, _beta, _gamma) *
      ElasticityTensorTools::elasticJacobian(_elasticity[_qp], _component, _component,
                                            _grad_test[_i][_qp], _grad_phi[_j][_qp]));
}

Real
CDPVelocityStressDivergence::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_zeta[_qp] == 0.0)
    return DynamicStressDivergenceTensors::computeQpOffDiagJacobian(jvar);
  if (_dt <= 0)
    return 0;
  for (unsigned int c = 0; c < _ndisp; ++c)
    if (jvar == _disp_var[c])
      return (1 + _alpha) * (StressDivergenceTensors::computeQpOffDiagJacobian(jvar) +
          _zeta[_qp] * ElasticVelocityRate::displacementDerivative(_dt, _beta, _gamma) *
          ElasticityTensorTools::elasticJacobian(_elasticity[_qp], _component, c,
                                                _grad_test[_i][_qp], _grad_phi[_j][_qp]));
  return (1 + _alpha) * StressDivergenceTensors::computeQpOffDiagJacobian(jvar);
}
