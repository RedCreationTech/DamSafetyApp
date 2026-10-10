#pragma once
#include "DynamicStressDivergenceTensors.h"

class CDPVelocityStressDivergence : public DynamicStressDivergenceTensors
{
public:
  static InputParameters validParams();
  CDPVelocityStressDivergence(const InputParameters & parameters);
protected:
  Real computeQpResidual() override;
  Real computeQpJacobian() override;
  Real computeQpOffDiagJacobian(unsigned int jvar) override;
  const Real _beta, _gamma;
  const MaterialProperty<RankFourTensor> & _elasticity;
  const MaterialProperty<RankTwoTensor> & _damping_stress;
};
