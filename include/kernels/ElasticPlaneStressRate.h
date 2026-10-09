#pragma once
#include "Kernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

class ElasticPlaneStressRate : public Kernel
{
public:
  static InputParameters validParams();
  ElasticPlaneStressRate(const InputParameters & parameters);
protected:
  Real computeQpResidual() override;
  Real computeQpJacobian() override;
  Real computeQpOffDiagJacobian(unsigned int jvar) override;
  const Real _beta, _gamma;
  const MaterialProperty<RankFourTensor> & _elasticity;
  const MaterialProperty<RankTwoTensor> & _rate_stress;
  std::vector<unsigned int> _displacements;
};
