#include "AcceptedIntervalStressDivergence.h"
#include "StressHistoryInterval.h"
#include <iomanip>

registerMooseObject("DamSafetyApp", AcceptedIntervalStressDivergence);

InputParameters
AcceptedIntervalStressDivergence::validParams()
{
  auto params = DynamicStressDivergenceTensors::validParams();
  params.addClassDescription("HHT stress damping with the previous stress rate divided by its previous accepted interval. Current-stress tangent and all other terms remain upstream.");
  params.addParam<bool>("history_diagnostics", false, "Log trial time/dt/previous accepted dt for small-model validation only");
  return params;
}

AcceptedIntervalStressDivergence::AcceptedIntervalStressDivergence(const InputParameters & parameters)
  : DynamicStressDivergenceTensors(parameters),
    _history_diagnostics(getParam<bool>("history_diagnostics"))
{
}

void
AcceptedIntervalStressDivergence::timestepSetup()
{
  DynamicStressDivergenceTensors::timestepSetup();
  if (_history_diagnostics)
    _console << std::setprecision(17) << "STRESS_HISTORY_INTERVAL time=" << _t
             << " dt=" << _dt << " previous_dt=" << _dt_old << " step=" << _t_step
             << std::endl;
}

Real
AcceptedIntervalStressDivergence::computeQpResidual()
{
  const Real original = DynamicStressDivergenceTensors::computeQpResidual();
  if (_static_initialization && _t == _dt)
    return original;
  const Real factor = StressHistoryInterval::correction(_alpha, _zeta[_qp], _dt, _dt_old, _t_step);
  if (factor == 0)
    return original;
  const auto history = _stress_older[_qp] - _stress_old[_qp];
  Real extra = history.row(_component) * _grad_test[_i][_qp];
  if (_volumetric_locking_correction)
    extra += history.trace() / 3.0 *
             (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component));
  // old/older material states are fixed during Newton: inherited Jacobians remain valid.
  return original + factor * extra;
}
