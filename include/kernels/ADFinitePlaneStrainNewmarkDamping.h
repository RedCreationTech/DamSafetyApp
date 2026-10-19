#pragma once

#include "ADKernel.h"
#include "RankFourTensor.h"

class ADFinitePlaneStrainNewmarkDamping : public ADKernel
{
public:
  static InputParameters validParams();
  ADFinitePlaneStrainNewmarkDamping(const InputParameters & parameters);
  void initialSetup() override;

protected:
  ADReal computeQpResidual() override;
  const unsigned int _component;
  const Real _beta, _gamma, _alpha, _youngs, _poisson, _expected_zeta;
  const MaterialProperty<Real> & _zeta;
  const MaterialProperty<RankFourTensor> & _elasticity;
  std::vector<const ADVariableGradient *> _gradients;
  std::vector<const VariableGradient *> _old_gradients, _old_velocities, _old_accelerations;
};
