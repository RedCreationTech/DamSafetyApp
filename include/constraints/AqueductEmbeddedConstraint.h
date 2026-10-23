#pragma once
#include "EqualValueEmbeddedConstraint.h"

// Use the audited source host map and one owning truss block per node.
class AqueductEmbeddedConstraint : public EqualValueEmbeddedConstraint
{
public:
  static InputParameters validParams();
  AqueductEmbeddedConstraint(const InputParameters &);
protected:
  void prepareSecondaryToPrimaryMap() override;
};
