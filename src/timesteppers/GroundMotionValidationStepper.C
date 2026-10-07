#include "GroundMotionValidationStepper.h"
#include <iomanip>
registerMooseObject("DamSafetyApp", GroundMotionValidationStepper);
InputParameters GroundMotionValidationStepper::validParams()
{
  auto params = TimeSequenceStepper::validParams();
  params.addClassDescription("Small-model validation: one deliberately rejected converged trial and full-precision accepted times. Never use to run a production earthquake.");
  params.addRequiredParam<bool>("validation_only", "Must be true to explicitly enable this diagnostic");
  params.addParam<Real>("reject_time", -1, "Reject first converged trial at or after this time; negative disables");
  return params;
}
GroundMotionValidationStepper::GroundMotionValidationStepper(const InputParameters & parameters)
  : TimeSequenceStepper(parameters), _reject_time(getParam<Real>("reject_time"))
{
  if (!getParam<bool>("validation_only")) paramError("validation_only", "Diagnostic requires explicit opt-in");
}
void GroundMotionValidationStepper::step()
{
  TimeStepper::step();
  if (_converged && !_rejected && _reject_time >= 0 && _time >= _reject_time)
  {
    _rejected = true;
    _converged = false;
    ++_failure_count;
    _console << std::setprecision(17) << "GROUND_MOTION_REJECT time=" << _time << " dt=" << _dt << '\n';
  }
}
void GroundMotionValidationStepper::acceptStep()
{
  TimeSequenceStepper::acceptStep();
  _console << std::setprecision(17) << "GROUND_MOTION_ACCEPT time=" << _time << " dt=" << _dt << '\n';
}
