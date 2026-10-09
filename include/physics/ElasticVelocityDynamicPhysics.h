#pragma once
#include "DynamicSolidMechanicsPhysics.h"

class ElasticVelocityDynamicPhysics : public DynamicSolidMechanicsPhysics
{
public:
  static InputParameters validParams();
  ElasticVelocityDynamicPhysics(const InputParameters & parameters);
  void act() override;
protected:
  std::string getKernelType() override;
  InputParameters getKernelParameters(std::string type) override;
  const bool _enabled;
};
