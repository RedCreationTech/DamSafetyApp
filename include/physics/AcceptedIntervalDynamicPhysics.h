#pragma once
#include "DynamicSolidMechanicsPhysics.h"

class AcceptedIntervalDynamicPhysics : public DynamicSolidMechanicsPhysics
{
public:
  static InputParameters validParams() { return DynamicSolidMechanicsPhysics::validParams(); }
  AcceptedIntervalDynamicPhysics(const InputParameters & parameters) : DynamicSolidMechanicsPhysics(parameters) {}
protected:
  std::string getKernelType() override
  {
    if (DynamicSolidMechanicsPhysics::getKernelType() != "DynamicStressDivergenceTensors")
      mooseError("AcceptedIntervalDynamicPhysics requires Cartesian dynamic mechanics");
    return "AcceptedIntervalStressDivergence";
  }
};
