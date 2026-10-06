#pragma once
#include "StressDivergenceTensors.h"
#include "UniformHex8.h"
#include <unordered_map>

class UniformHex8StressDivergence : public StressDivergenceTensors
{
public:
  static InputParameters validParams();
  UniformHex8StressDivergence(const InputParameters & parameters);
protected:
  Real computeQpResidual() override;
  Real computeQpJacobian() override;
  Real computeQpOffDiagJacobian(unsigned int jvar) override;
  const UniformHex8::Operators & operators();
  Real stiffness(unsigned int column);
  const Real _initial_shear_modulus;
  const Real _hourglass_factor;
  MooseVariable * _nodal_displacement;
  std::unordered_map<dof_id_type, UniformHex8::Operators> _operators;
};
