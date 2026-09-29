#pragma once

#include "ADKernel.h"
#include "Quad4Hourglass.h"

class ADQuad4HourglassControl : public ADKernel
{
public:
  static InputParameters validParams();
  ADQuad4HourglassControl(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  const Real _coefficient;
  const Real _shear_modulus;
  const ADDofValues & _dof_values;
};
