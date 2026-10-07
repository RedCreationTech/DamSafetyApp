#pragma once
#include "DynamicStressDivergenceTensors.h"

class AcceptedIntervalStressDivergence : public DynamicStressDivergenceTensors
{
public:
  static InputParameters validParams();
  AcceptedIntervalStressDivergence(const InputParameters & parameters);
  void timestepSetup() override;
protected:
  Real computeQpResidual() override;
  const bool _history_diagnostics;
};
