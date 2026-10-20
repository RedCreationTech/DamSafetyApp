#include "CDPObservationTimeStepper.h"
#include "CDPObservationProblem.h"
#include "NonlinearSystem.h"
#include "libmesh/petsc_nonlinear_solver.h"
registerMooseObject("DamSafetyApp", CDPObservationTimeStepper);
InputParameters CDPObservationTimeStepper::validParams() { return TimeSequenceStepper::validParams(); }
void CDPObservationTimeStepper::step()
{
  auto * p=dynamic_cast<CDPObservationProblem *>(&_fe_problem);
  if(p) p->attemptEvent("attempt_enter");
  TimeSequenceStepper::step();
  if(p) p->attemptEvent(converged()?"solve_converged":"solve_failed");
}
void CDPObservationTimeStepper::rejectStep()
{
  auto * p=dynamic_cast<CDPObservationProblem *>(&_fe_problem);
  if(p) p->attemptEvent("reject_before");
  TimeSequenceStepper::rejectStep();
  if(p) p->attemptEvent("reject_after");
}
void CDPObservationTimeStepper::acceptStep()
{
  auto * p=dynamic_cast<CDPObservationProblem *>(&_fe_problem);
  if(p) p->attemptEvent("accepted_end");
  TimeSequenceStepper::acceptStep();
}
