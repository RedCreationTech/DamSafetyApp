#pragma once
#include "ADMortarConstraint.h"
#include "SmallSlidingCohesiveContact.h"
class SmallSlidingHardContact : public ADMortarConstraint
{
public:
  static InputParameters validParams();
  SmallSlidingHardContact(const InputParameters & p);
  void post() override;
protected:
  ADReal computeQpResidual(Moose::MortarType) override { return 0; }
  void computeResidual(Moose::MortarType) override {}
  void computeJacobian(Moose::MortarType) override {}
  const Real _c;
  const SmallSlidingCohesiveContact & _contact;
};
