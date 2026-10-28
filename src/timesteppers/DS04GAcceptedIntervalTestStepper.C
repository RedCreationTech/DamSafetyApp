#include "DS04GAcceptedIntervalTestStepper.h"
#include "FEProblemBase.h"
#include "MooseMesh.h"
#include <iomanip>

registerMooseObject("DamSafetyApp", DS04GAcceptedIntervalTestStepper);
InputParameters DS04GAcceptedIntervalTestStepper::validParams()
{
  auto p=ConstantDT::validParams();
  p.addClassDescription("TEST ONLY: one injected failure to verify real accepted-dt semantics");
  return p;
}
DS04GAcceptedIntervalTestStepper::DS04GAcceptedIntervalTestStepper(const InputParameters & p)
  : ConstantDT(p), _trace("accepted_interval_executioner.csv")
{
  _trace << std::setprecision(17) << "event,step,time,h,h_previous,last_accepted_h\n";
}
void DS04GAcceptedIntervalTestStepper::record(const char * event)
{
  _trace << event << ',' << _t_step << ',' << _time << ',' << _dt << ','
         << _fe_problem.dtOld() << ',' << _last_accepted << '\n';
  if (!_trace) mooseError("Could not write accepted interval trace");
}
void DS04GAcceptedIntervalTestStepper::step()
{
  if (comm().size()!=1 || _end_time>.01 || _accepted>=8 ||
      _fe_problem.mesh().getMesh().n_elem()!=1)
    mooseError("Accepted interval tester is bounded to serial one element and at most .01s/8 steps");
  record("attempt");
  if (_accepted && _fe_problem.dtOld()!=_last_accepted)
    mooseError("Previous accepted dt was changed by a rejected attempt");
  if (_cutback_accepted && _dt>_last_accepted) _after_cutback=true;
  TimeStepper::step();
  if (_t_step==2 && !_injected)
  {
    if (!_converged) mooseError("The bounded test solve unexpectedly failed before injection");
    _injected=true; _converged=false; ++_failure_count;
    record("injected_failure");
  }
}
void DS04GAcceptedIntervalTestStepper::acceptStep()
{
  record("accept");
  if (_injected && _dt<_last_accepted) _cutback_accepted=true;
  _last_accepted=_dt; ++_accepted;
  ConstantDT::acceptStep();
}
void DS04GAcceptedIntervalTestStepper::rejectStep()
{
  const Real hp=_fe_problem.dtOld();
  ConstantDT::rejectStep();
  if (_fe_problem.dtOld()!=hp) mooseError("rejectStep modified accepted dt");
  record("reject_restore");
}
void DS04GAcceptedIntervalTestStepper::postExecute()
{
  ConstantDT::postExecute();
  if (!_injected || !_cutback_accepted || !_after_cutback || _accepted>8 || !_accepted)
    mooseError("Incomplete or unbounded accepted-interval test");
  mooseInfo("DS04G accepted interval executioner check complete");
}
