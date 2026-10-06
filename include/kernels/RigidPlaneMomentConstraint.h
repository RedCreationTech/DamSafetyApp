#pragma once
#include "NodalScalarKernel.h"
#include <array>

// Free rotation equations: - integral (X-X_RP) x lambda dA = 0.
// moment_file stores exact reference-face shape integrals, not force results.
class RigidPlaneMomentConstraint : public NodalScalarKernel
{
public:
  static InputParameters validParams();
  RigidPlaneMomentConstraint(const InputParameters & parameters);
  void computeResidual() override;
  void computeJacobian() override;
protected:
  Real coefficient(unsigned int rotation, unsigned int traction, unsigned int node) const;
  const std::vector<unsigned int> _traction_numbers;
  const std::vector<const VariableValue *> _tractions;
  std::vector<RealVectorValue> _moments;
};
