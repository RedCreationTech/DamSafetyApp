#pragma once
#include "ADMortarLagrangeConstraint.h"
#include "SmallSlidingCohesiveContact.h"
class SmallSlidingContactTraction : public ADMortarLagrangeConstraint
{
public:
  static InputParameters validParams();
  SmallSlidingContactTraction(const InputParameters & p);
protected:
  ADReal computeQpResidual(Moose::MortarType type) override;
  const unsigned int _component;
  const SmallSlidingCohesiveContact & _contact;
};
