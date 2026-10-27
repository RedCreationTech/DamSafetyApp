#pragma once
#include "DynamicStressDivergenceTensors.h"

/** Opt-in Jacobian only. Residual delegates to the exact original static/dynamic kernel. */
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
  const bool _sp1_dynamic;
  const MaterialProperty<RankTwoTensor> & _sp1_F;
  const MaterialProperty<RankTwoTensor> & _sp1_Fold;
};
