#pragma once

#include "DynamicSolidMechanicsPhysics.h"

class FinitePlaneStrainNewmarkPhysics : public DynamicSolidMechanicsPhysics
{
public:
  static InputParameters validParams();
  FinitePlaneStrainNewmarkPhysics(const InputParameters & parameters);
  void act() override;

protected:
  InputParameters getKernelParameters(std::string type) override;
  const bool _enabled;
  const bool _sp1;
};
