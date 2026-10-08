#include "ADQuad4HourglassControl.h"
#include "StressHistoryInterval.h"
#include <iomanip>

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
  params.addParam<MaterialPropertyName>(
      "zeta", 0.0, "Stiffness-proportional Rayleigh damping coefficient");
  params.addParam<Real>("alpha", 0.0, "HHT time integration parameter");
  params.addParam<bool>("use_previous_accepted_interval", false,
                       "Use the previous accepted timestep for the previous hourglass rate; opt-in diagnostic, default retains original behavior");
  params.addParam<bool>("history_diagnostics", false,
                       "Log trial time/dt/previous accepted dt for small-model validation only");
  return params;
}

ADQuad4HourglassControl::ADQuad4HourglassControl(const InputParameters & parameters)
  : ADKernel(parameters),
    _coefficient(getParam<Real>("coefficient")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _zeta(getMaterialProperty<Real>("zeta")),
    _alpha(getParam<Real>("alpha")),
    _previous_accepted_interval(getParam<bool>("use_previous_accepted_interval")),
    _history_diagnostics(getParam<bool>("history_diagnostics")),
    _dof_values(_var.adDofValues()),
    _dof_values_old(_var.dofValuesOld()),
    _dof_values_older(_var.dofValuesOlder())
{
}

void
ADQuad4HourglassControl::timestepSetup()
{
  ADKernel::timestepSetup();
  if (_history_diagnostics)
    _console << std::setprecision(17) << "HOURGLASS_HISTORY_INTERVAL time=" << _t
             << " dt=" << _dt << " previous_dt=" << _dt_old << " step=" << _t_step
             << std::endl;
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
  Real amplitude_old = 0.0;
  Real amplitude_older = 0.0;
  for (unsigned int i = 0; i < 4; ++i)
  {
    amplitude += gamma[i] * _dof_values[i];
    amplitude_old += gamma[i] * _dof_values_old[i];
    amplitude_older += gamma[i] * _dof_values_older[i];
  }

  if (_dt > 0)
  {
    amplitude = amplitude * (1.0 + _alpha + (1.0 + _alpha) * _zeta[_qp] / _dt) -
                amplitude_old * (_alpha + (1.0 + 2.0 * _alpha) * _zeta[_qp] / _dt) +
                amplitude_older * _alpha * _zeta[_qp] / _dt;
    if (_previous_accepted_interval)
      // The accepted old/older modes are fixed in Newton, so the current AD tangent is unchanged.
      amplitude += StressHistoryInterval::correction(_alpha, _zeta[_qp], _dt, _dt_old, _t_step) *
                   (amplitude_older - amplitude_old);
  }

  const Real length_squared = Quad4Hourglass::characteristicLengthSquared(nodes);
  return _coefficient * _shear_modulus / length_squared * gamma[_i] * amplitude;
}
