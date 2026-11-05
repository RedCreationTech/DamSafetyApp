#include "CDPObservationConstantDT.h"

#include "CDPObservationProblem.h"

registerMooseObject("DamSafetyApp", CDPObservationConstantDT);

InputParameters
CDPObservationConstantDT::validParams()
{
  return ConstantDT::validParams();
}

CDPObservationConstantDT::CDPObservationConstantDT(const InputParameters & parameters)
  : ConstantDT(parameters)
{
}

void
CDPObservationConstantDT::step()
{
  auto * problem = dynamic_cast<CDPObservationProblem *>(&_fe_problem);
  const bool observe = problem && problem->observationEnabled();
  if (observe)
    problem->attemptEvent("attempt_enter");
  try
  {
    ConstantDT::step();
  }
  catch (...)
  {
    if (observe)
      problem->attemptEvent("solve_exception");
    throw;
  }
  if (observe)
    problem->attemptEvent(converged() ? "solve_converged" : "solve_failed");
}

void
CDPObservationConstantDT::rejectStep()
{
  auto * problem = dynamic_cast<CDPObservationProblem *>(&_fe_problem);
  const bool observe = problem && problem->observationEnabled();
  if (observe)
    problem->attemptEvent("reject_before");
  try
  {
    ConstantDT::rejectStep();
  }
  catch (...)
  {
    if (observe)
      problem->attemptEvent("reject_exception");
    throw;
  }
  if (observe)
    problem->attemptEvent("reject_after");
}

void
CDPObservationConstantDT::acceptStep()
{
  auto * problem = dynamic_cast<CDPObservationProblem *>(&_fe_problem);
  const bool observe = problem && problem->observationEnabled();
  if (observe)
    problem->attemptEvent("accepted_end");
  try
  {
    ConstantDT::acceptStep();
  }
  catch (...)
  {
    if (observe)
      problem->attemptEvent("accept_exception");
    throw;
  }
  if (observe)
    problem->attemptEvent("stepper_accept_after");
}
