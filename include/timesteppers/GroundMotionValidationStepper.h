#pragma once
#include "TimeSequenceStepper.h"
/** Opt-in diagnostic only: solve and then reject one trial to exercise rollback. */
class GroundMotionValidationStepper : public TimeSequenceStepper
{
public:
  static InputParameters validParams();
  GroundMotionValidationStepper(const InputParameters & parameters);
  void step() override;
  void acceptStep() override;
private:
  const Real _reject_time;
  bool _rejected = false;
};
