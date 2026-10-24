#pragma once
#include "Function.h"
#include "IntegratedGroundMotion.h"
#include <memory>

/** Independent U/V/A from one original acceleration record. */
class IntegratedAcceleration : public Function
{
public:
  static InputParameters validParams();
  IntegratedAcceleration(const InputParameters & parameters);
  Real value(Real t, const Point & p) const override;
  Real timeDerivative(Real t, const Point & p) const override;
private:
  const MooseEnum _component;
  const Real _start_time;
  const Real _initial_displacement;
  std::unique_ptr<IntegratedGroundMotion> _motion;
};
