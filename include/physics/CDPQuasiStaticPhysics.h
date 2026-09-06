#pragma once
#include "QuasiStaticSolidMechanicsPhysics.h"
class CDPQuasiStaticPhysics : public QuasiStaticSolidMechanicsPhysics
{
public:
  static InputParameters validParams() { return QuasiStaticSolidMechanicsPhysics::validParams(); }
  CDPQuasiStaticPhysics(const InputParameters & p) : QuasiStaticSolidMechanicsPhysics(p) {}
protected:
  std::string getKernelType() override
  {
    const auto original=QuasiStaticSolidMechanicsPhysics::getKernelType();
    if (original!="StressDivergenceTensors")
      mooseError("CDPQuasiStaticPhysics supports Cartesian small-strain mechanics only");
    return "CDPStressDivergenceTensors";
  }
};
