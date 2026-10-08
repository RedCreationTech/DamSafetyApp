#pragma once

#include "ADKernel.h"
#include "Quad4Hourglass.h"

class ADQuad4HourglassControl : public ADKernel
{
public:
  static InputParameters validParams();
  ADQuad4HourglassControl(const InputParameters & parameters);

protected:
  void timestepSetup() override;
  ADReal computeQpResidual() override;

  const Real _coefficient;
  const Real _shear_modulus;
  const MaterialProperty<Real> & _zeta;
  const Real _alpha;
  const bool _previous_accepted_interval;
  const bool _history_diagnostics;
  const MooseArray<ADReal> & _dof_values;
  const VariableValue & _dof_values_old;
  const VariableValue & _dof_values_older;
};
