#include "ADQuad4ConsistentMassCorrection.h"

#include "libmesh/elem.h"
#include "Function.h"

registerMooseObject("DamSafetyApp", ADQuad4ConsistentMassCorrection);

InputParameters
ADQuad4ConsistentMassCorrection::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription(
      "Adds the exact QUAD4 consistent-mass contribution omitted by one-point quadrature.");
  params.addRequiredCoupledVar("velocity", "Newmark velocity variable");
  params.addRequiredCoupledVar("acceleration", "Newmark acceleration variable");
  params.addRequiredRangeCheckedParam<Real>("beta", "beta > 0", "Newmark beta");
  params.addRequiredRangeCheckedParam<Real>("gamma", "gamma > 0", "Newmark gamma");
  params.addParam<Real>("alpha", 0.0, "HHT alpha");
  params.addParam<MaterialPropertyName>("density", "density", "Density material property");
  params.addParam<MaterialPropertyName>(
      "eta", 0.0, "Mass-proportional Rayleigh damping material property");
  params.addParam<FunctionName>("ground_acceleration", "Optional relative-coordinate ground acceleration");
  params.addParam<FunctionName>("ground_velocity", "Optional relative-coordinate ground velocity");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

ADQuad4ConsistentMassCorrection::ADQuad4ConsistentMassCorrection(
    const InputParameters & parameters)
  : ADKernel(parameters),
    _ground_acceleration(isParamValid("ground_acceleration") ? &getFunction("ground_acceleration") : nullptr),
    _ground_velocity(isParamValid("ground_velocity") ? &getFunction("ground_velocity") : nullptr),
    _beta(getParam<Real>("beta")),
    _gamma(getParam<Real>("gamma")),
    _alpha(getParam<Real>("alpha")),
    _density(getMaterialProperty<Real>("density")),
    _eta(getMaterialProperty<Real>("eta")),
    _dof_values(_var.adDofValues()),
    _dof_values_old(_var.dofValuesOld()),
    _velocity_old(coupledDofValuesOld("velocity")),
    _acceleration_old(coupledDofValuesOld("acceleration"))
{
}

ADReal
ADQuad4ConsistentMassCorrection::computeQpResidual()
{
  if (_qrule->n_points() != 1)
    mooseError(name(), " requires one-point quadrature");
  if (_current_elem->type() != QUAD4 || _dof_values.size() != 4)
    mooseError(name(), " supports first-order QUAD4 elements only");
  if (_dt <= 0)
    return 0.0;

  Quad4Hourglass::Nodes nodes;
  for (unsigned int i = 0; i < 4; ++i)
  {
    nodes[i][0] = _current_elem->point(i)(0);
    nodes[i][1] = _current_elem->point(i)(1);
  }
  const auto correction = Quad4Hourglass::consistentMassCorrection(nodes);

  // The correction row sum is generally nonzero on distorted quads. Restore
  // this part of M*1*ground_load as well as the relative-coordinate mass.
  if (bool(_ground_acceleration) != bool(_ground_velocity))
    mooseError(name(), " requires both ground_acceleration and ground_velocity, or neither");
  const Real ground_load = _ground_acceleration ?
      _ground_acceleration->value(_t, _q_point[_qp]) + _eta[_qp] *
      ((1.0 + _alpha) * _ground_velocity->value(_t, _q_point[_qp]) -
       _alpha * _ground_velocity->value(_t - _dt, _q_point[_qp])) : 0.0;
  ADReal force = 0.0;
  for (unsigned int j = 0; j < 4; ++j)
  {
    const ADReal acceleration =
        ((_dof_values[j] - _dof_values_old[j]) / (_dt * _dt) - _velocity_old[j] / _dt -
         _acceleration_old[j] * (0.5 - _beta)) /
        _beta;
    const ADReal velocity = _velocity_old[j] + _dt * (1.0 - _gamma) * _acceleration_old[j] +
                            _gamma * _dt * acceleration;
    force += correction[_i][j] *
             (acceleration + ground_load +
              _eta[_qp] * ((1.0 + _alpha) * velocity - _alpha * _velocity_old[j]));
  }
  return _density[_qp] * force / Quad4Hourglass::area(nodes);
}
