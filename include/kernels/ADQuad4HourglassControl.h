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
  const MaterialProperty<Real> & _zeta;
  const Real _alpha;
  const MooseArray<ADReal> & _dof_values;
  const DofValues & _dof_values_old;
  const DofValues & _dof_values_older;
};
