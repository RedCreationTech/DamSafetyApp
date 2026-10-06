#pragma once
#include "AuxKernel.h"
#include "SmallSlidingCohesiveContact.h"
class SmallSlidingContactAux : public AuxKernel
{
public:
  static InputParameters validParams();
  SmallSlidingContactAux(const InputParameters & p);
protected:
  Real computeValue() override;
  const SmallSlidingCohesiveContact & _contact;
  const unsigned int _quantity;
};
