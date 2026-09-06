#pragma once
#include "StressDivergenceTensors.h"
class CDPStressDivergenceTensors : public StressDivergenceTensors
{
public:
  static InputParameters validParams();
  CDPStressDivergenceTensors(const InputParameters & p);
protected:
  Real computeQpJacobian() override;
  Real computeQpOffDiagJacobian(unsigned int jvar) override;
  Real bbarJacobian(unsigned int component) const;
};
