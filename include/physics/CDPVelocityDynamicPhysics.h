#pragma once
#include "DynamicSolidMechanicsPhysics.h"

class CDPVelocityDynamicPhysics : public DynamicSolidMechanicsPhysics
{
public:
  static InputParameters validParams();
  CDPVelocityDynamicPhysics(const InputParameters & parameters);
  void act() override;
protected:
  std::string getKernelType() override;
  InputParameters getKernelParameters(std::string type) override;
  const bool _enabled, _plane_stress_projection, _complete_coupling;
};
