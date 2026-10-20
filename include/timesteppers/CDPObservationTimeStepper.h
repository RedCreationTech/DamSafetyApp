#pragma once
#include "TimeSequenceStepper.h"
/** Delegate unchanged sequence/cutback behavior, reading state at attempt boundaries. */
class CDPObservationTimeStepper : public TimeSequenceStepper
{
public:
  static InputParameters validParams();
  CDPObservationTimeStepper(const InputParameters & p) : TimeSequenceStepper(p) {}
  void step() override;
  void rejectStep() override;
  void acceptStep() override;
};
