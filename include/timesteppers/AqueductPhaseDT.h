#pragma once
#include "FunctionDT.h"

/** Source static/dynamic timestep bounds with FunctionDT phase knots. */
class AqueductPhaseDT : public FunctionDT
{
public:
  static InputParameters validParams();
  AqueductPhaseDT(const InputParameters &);
  bool constrainStep(Real &) override;
  void acceptStep() override;
protected:
  Real computeDT() override;
  Real computeFailedDT() override;
  const Real _static_end,_static_min,_static_max,_dynamic_min,_dynamic_max;
  const unsigned int _dynamic_retry_hold_steps;
  unsigned int & _hold_remaining;
  bool & _recovering_from_failure;
  Real & _hold_ceiling;
};
