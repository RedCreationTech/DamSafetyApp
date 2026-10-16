#pragma once
#include "FEProblem.h"

/** Opt-in static-to-dynamic kinematic initialization; never edits material state. */
class StaticRestartHistoryProblem : public FEProblem
{
public:
  static InputParameters validParams();
  StaticRestartHistoryProblem(const InputParameters & parameters);
  void initialSetup() override;

private:
  bool & _static_history_initialized;
};
