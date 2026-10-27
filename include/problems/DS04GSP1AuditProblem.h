#pragma once
#include "FEProblem.h"

/** Test-only frozen-state assembly audit; admitted only for one serial QUAD4. */
class DS04GSP1AuditProblem : public FEProblem
{
public:
  static InputParameters validParams();
  DS04GSP1AuditProblem(const InputParameters &);
  void initialSetup() override;
};
