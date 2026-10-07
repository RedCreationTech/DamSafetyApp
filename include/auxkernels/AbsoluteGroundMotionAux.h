#pragma once
#include "AuxKernel.h"
class AbsoluteGroundMotionAux : public AuxKernel
{
public:
  static InputParameters validParams();
  AbsoluteGroundMotionAux(const InputParameters & parameters);
protected:
  Real computeValue() override;
  const VariableValue & _relative;
  const Function & _ground;
};
