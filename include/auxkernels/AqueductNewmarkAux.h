#pragma once
#include "NewmarkAccelAux.h"
#include "NewmarkVelAux.h"

/** Physical motion rates exclude deformation accumulated during quasi-static loading. */
class AqueductNewmarkAcceleration : public NewmarkAccelAux
{
public:
  static InputParameters validParams();
  AqueductNewmarkAcceleration(const InputParameters &);
protected:
  Real computeValue() override;
  const Real _static_end;
};

class AqueductNewmarkVelocity : public NewmarkVelAux
{
public:
  static InputParameters validParams();
  AqueductNewmarkVelocity(const InputParameters &);
protected:
  Real computeValue() override;
  const Real _static_end;
};
