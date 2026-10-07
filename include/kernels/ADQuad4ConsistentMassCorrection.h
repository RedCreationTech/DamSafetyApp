#pragma once

#include "ADKernel.h"
#include "Quad4Hourglass.h"

class ADQuad4ConsistentMassCorrection : public ADKernel
{
public:
  static InputParameters validParams();
  ADQuad4ConsistentMassCorrection(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  const Function * const _ground_acceleration;
  const Function * const _ground_velocity;
  const Real _beta;
  const Real _gamma;
  const Real _alpha;
  const MaterialProperty<Real> & _density;
  const MaterialProperty<Real> & _eta;
  const MooseArray<ADReal> & _dof_values;
  const VariableValue & _dof_values_old;
  const VariableValue & _velocity_old;
  const VariableValue & _acceleration_old;
};
