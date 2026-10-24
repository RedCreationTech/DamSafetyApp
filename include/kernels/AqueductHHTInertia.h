#pragma once
#include "InertialForce.h"

/** Consistent solid/truss inertia for normalized whole-system HHT. */
class AqueductHHTInertia : public InertialForce
{
public:
  static InputParameters validParams();
  AqueductHHTInertia(const InputParameters &);
protected:
  Real computeQpResidual() override;
  Real computeQpJacobian() override;
  Real _hht_mass_scale,_hht_start_time;
};
