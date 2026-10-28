#pragma once
#include "FEProblem.h"

/** TEST ONLY: one serial QUAD4, prescribed histories, no accepted time step. */
class DS04GSP1DynamicAuditProblem : public FEProblem
{
public:
  static InputParameters validParams();
  DS04GSP1DynamicAuditProblem(const InputParameters &);
  void initialSetup() override;
};
