#pragma once
#include "FunctionDT.h"

/** Source static/dynamic timestep bounds with FunctionDT phase knots. */
class AqueductPhaseDT : public FunctionDT
{
public:
  static InputParameters validParams();
  AqueductPhaseDT(const InputParameters &);
  bool constrainStep(Real &) override;
protected:
  Real computeFailedDT() override;
  const Real _static_end,_static_min,_static_max,_dynamic_min,_dynamic_max;
};
