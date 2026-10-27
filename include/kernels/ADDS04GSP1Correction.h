#pragma once
#include "ADKernel.h"
#include "DS04GSP1.h"
#include "RankFourTensor.h"

class ADDS04GSP1Correction : public ADKernel
{
public:
  static InputParameters validParams();
  ADDS04GSP1Correction(const InputParameters &);
  void initialSetup() override;

protected:
  ADReal computeQpResidual() override;
  const MooseEnum _term;
  const unsigned _component;
  const Real _youngs, _nu, _rho, _eta, _zeta, _alpha, _beta, _gamma, _static_end, _gravity;
  const MaterialProperty<Real> &_density;
  const MaterialProperty<RankFourTensor> &_elasticity;
  std::array<const ADVariableValue *, 2> _u;
  const VariableValue &_old;
  const VariableValue &_older;
  const VariableValue &_v_old;
  const VariableValue &_a_old;
};
