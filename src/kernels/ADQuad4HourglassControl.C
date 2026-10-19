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
  params.addParam<MooseEnum>(
      "elastic_operator", MooseEnum("legacy gradient_raw_conjugate", "legacy"),
      "Opt-in energy-checked elastic geometry/amplitude; retains the original damping operator");
  params.addParam<bool>("use_previous_accepted_interval", false,
                       "Use the previous accepted timestep for the previous hourglass rate; opt-in diagnostic, default retains original behavior");
  params.addParam<bool>("history_diagnostics", false,
                       "Log trial time/dt/previous accepted dt for small-model validation only");
  params.addParam<MooseEnum>("damping_rate", MooseEnum("legacy_secant newmark_velocity", "legacy_secant"),
      "Opt-in Newmark velocity rate for the separate hourglass damping; no spatial coefficient change");
  params.addCoupledVar("velocity", "Nodal Newmark velocity; required only by newmark_velocity");
  params.addCoupledVar("acceleration", "Nodal Newmark acceleration; required only by newmark_velocity");
  params.addRangeCheckedParam<Real>("newmark_beta", 0.25, "newmark_beta > 0", "Newmark beta");
  params.addRangeCheckedParam<Real>("newmark_gamma", 0.5, "newmark_gamma > 0", "Newmark gamma");
  return params;
}

ADQuad4HourglassControl::ADQuad4HourglassControl(const InputParameters & parameters)
  : ADKernel(parameters),
    _coefficient(getParam<Real>("coefficient")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _zeta(getMaterialProperty<Real>("zeta")),
    _alpha(getParam<Real>("alpha")),
    _gradient_raw_conjugate(getParam<MooseEnum>("elastic_operator") == "gradient_raw_conjugate"),
    _previous_accepted_interval(getParam<bool>("use_previous_accepted_interval")),
    _history_diagnostics(getParam<bool>("history_diagnostics")),
    _newmark_velocity_rate(getParam<MooseEnum>("damping_rate") == "newmark_velocity"),
    _newmark_beta(getParam<Real>("newmark_beta")),
    _newmark_gamma(getParam<Real>("newmark_gamma")),
    _velocity_old(_newmark_velocity_rate && isCoupled("velocity")
                      ? &coupledDofValuesOld("velocity") : nullptr),
    _acceleration_old(_newmark_velocity_rate && isCoupled("acceleration")
                      ? &coupledDofValuesOld("acceleration") : nullptr),
    _dof_values(_var.adDofValues()),
    _dof_values_old(_var.dofValuesOld()),
    _dof_values_older(_var.dofValuesOlder())
{
  if (_gradient_raw_conjugate && getParam<bool>("use_displaced_mesh"))
    paramError("use_displaced_mesh", "gradient_raw_conjugate requires reference geometry");
  if (_newmark_velocity_rate && (!_velocity_old || !_acceleration_old))
    paramError("damping_rate", "newmark_velocity requires coupled velocity and acceleration");
  if (_newmark_velocity_rate && getParam<bool>("use_displaced_mesh"))
    paramError("use_displaced_mesh", "newmark_velocity requires reference geometry");
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
  if (_newmark_velocity_rate)
  {
    ADReal rate = 0.0;
    if (_dt > 0)
    {
      if (_velocity_old->size() != 4 || _acceleration_old->size() != 4)
        mooseError(name(), " requires first-order nodal velocity and acceleration");
      for (unsigned int j = 0; j < 4; ++j)
        rate += gamma[j] * Quad4Hourglass::newmarkWeightedRate(
            _dof_values[j], _dof_values_old[j], (*_velocity_old)[j],
            (*_acceleration_old)[j], _dt, _newmark_beta, _newmark_gamma, _alpha);
    }
    const ADReal elastic = _dt > 0
        ? (1.0 + _alpha) * amplitude - _alpha * amplitude_old : amplitude;
    const Real elastic_stiffness = _gradient_raw_conjugate
        ? Quad4Hourglass::gradientRawStiffnessPerShear(nodes)
        : 1.0 / Quad4Hourglass::characteristicLengthSquared(nodes);
    // Keep the legacy spatial damping coefficient even when the elastic K is opt-in.
    return _coefficient * _shear_modulus * gamma[_i] *
           (elastic_stiffness * elastic + _zeta[_qp] * rate /
                Quad4Hourglass::characteristicLengthSquared(nodes));
  }
  if (_gradient_raw_conjugate)
  {
    const Real history_correction = _dt > 0 && _previous_accepted_interval
        ? StressHistoryInterval::correction(_alpha, _zeta[_qp], _dt, _dt_old, _t_step)
        : 0.0;
    return _coefficient * _shear_modulus * gamma[_i] *
           Quad4Hourglass::gradientRawResidualPerShear(
               nodes, _alpha, _dt, _zeta[_qp], history_correction,
               amplitude, amplitude_old, amplitude_older);
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
