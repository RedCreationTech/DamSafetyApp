#pragma once
#include "DynamicStressDivergenceTensors.h"

class ElasticVelocityStressDivergence : public DynamicStressDivergenceTensors
{
public:
  static InputParameters validParams();
  ElasticVelocityStressDivergence(const InputParameters & parameters);
protected:
  Real computeQpResidual() override;
  Real computeQpJacobian() override;
  Real computeQpOffDiagJacobian(unsigned int jvar) override;
  const Real _beta, _gamma;
  const unsigned int _rate_var;
  const MaterialProperty<RankFourTensor> & _elasticity;
  const MaterialProperty<RankTwoTensor> & _damping_stress;
};
