#pragma once

#include "ConstantDT.h"

class CDPObservationConstantDT : public ConstantDT
{
public:
  static InputParameters validParams();
  CDPObservationConstantDT(const InputParameters & parameters);

  void step() override;
  void rejectStep() override;
  void acceptStep() override;
};
