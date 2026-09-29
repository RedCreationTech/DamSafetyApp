#include "ADQuad4HourglassControl.h"

#include "libmesh/elem.h"

registerMooseObject("DamSafetyApp", ADQuad4HourglassControl);

InputParameters
ADQuad4HourglassControl::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription(
      "Adds an elastic stiffness penalty to the non-affine mode of a one-point QUAD4.");
  params.addRequiredRangeCheckedParam<Real>(
      "coefficient", "coefficient > 0", "Dimensionless hourglass stiffness coefficient");
  params.addRequiredRangeCheckedParam<Real>(
      "shear_modulus", "shear_modulus > 0", "Undamaged elastic shear modulus");
  return params;
}

ADQuad4HourglassControl::ADQuad4HourglassControl(const InputParameters & parameters)
  : ADKernel(parameters),
    _coefficient(getParam<Real>("coefficient")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _dof_values(_var.adDofValues())
{
}

ADReal
ADQuad4HourglassControl::computeQpResidual()
{
  if (_qrule->n_points() != 1)
    mooseError(name(), " requires one-point quadrature");
  if (_current_elem->type() != QUAD4 || _dof_values.size() != 4)
    mooseError(name(), " supports first-order QUAD4 elements only");

  Quad4Hourglass::Nodes nodes;
  for (unsigned int i = 0; i < 4; ++i)
  {
    nodes[i][0] = _current_elem->point(i)(0);
    nodes[i][1] = _current_elem->point(i)(1);
  }

  const auto gamma = Quad4Hourglass::mode(nodes);
  ADReal amplitude = 0.0;
  for (unsigned int i = 0; i < 4; ++i)
    amplitude += gamma[i] * _dof_values[i];

  const Real length_squared = Quad4Hourglass::characteristicLengthSquared(nodes);
  return _coefficient * _shear_modulus / length_squared * gamma[_i] * amplitude;
}
