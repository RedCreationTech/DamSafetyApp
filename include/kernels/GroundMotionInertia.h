#pragma once
#include "Kernel.h"
#include "FunctionInterface.h"

/** Ground-motion part of rho*(q_accel+Ag+eta*HHT(q_vel+Vg)). */
class GroundMotionInertia : public Kernel
{
public:
  static InputParameters validParams();
  GroundMotionInertia(const InputParameters & parameters);
protected:
  Real computeQpResidual() override;
  const Function & _acceleration;
  const Function & _velocity;
  const MaterialProperty<Real> & _density;
  const MaterialProperty<Real> & _eta;
  const Real _alpha;
};
