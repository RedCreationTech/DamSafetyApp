#pragma once
#include "DynamicStressDivergenceTensors.h"
#include "DS04GSP1.h"

/** Opt-in Jacobian only. Original residual by default; optional two-interval dynamic history. */
class DS04GSP1StressDivergence : public DynamicStressDivergenceTensors
{
public:
  static InputParameters validParams();
  DS04GSP1StressDivergence(const InputParameters &);

protected:
  void initialSetup() override;
  Real computeQpResidual() override;
  Real computeQpJacobian() override;
  Real computeQpOffDiagJacobian(unsigned) override;
  Real tangent(unsigned direction);
  DS04GSP1::Weights historyWeights() const;
  const bool _sp1_dynamic;
  const bool _two_interval;
  const MaterialProperty<RankTwoTensor> & _sp1_F;
  const MaterialProperty<RankTwoTensor> & _sp1_Fold;
};
