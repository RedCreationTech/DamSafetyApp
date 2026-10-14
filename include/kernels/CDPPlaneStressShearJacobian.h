#pragma once
#include "WeakPlaneStress.h"

// Adds the missing cross-block to the existing static WeakPlaneStress row.
class CDPPlaneStressShearJacobian : public WeakPlaneStress
{
public:
  static InputParameters validParams();
  CDPPlaneStressShearJacobian(const InputParameters & parameters);
protected:
  Real computeQpResidual() override { return 0.; }
  Real computeQpJacobian() override { return 0.; }
  Real computeQpOffDiagJacobian(unsigned int jvar) override;
  const bool _complete_coupling;
};
