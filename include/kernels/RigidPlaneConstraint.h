#pragma once
#include "Kernel.h"
#include "Function.h"

// Lagrange-multiplier equation u = U_RP(t) + theta x (X-X_RP).
// Applied on a lower-dimensional mesh sharing the original top-surface nodes.
class RigidPlaneConstraint : public Kernel
{
public:
  static InputParameters validParams();
  RigidPlaneConstraint(const InputParameters & parameters);
protected:
  Real computeQpResidual() override;
  Real computeQpJacobian() override { return 0.0; }
  Real computeQpOffDiagJacobian(unsigned int jvar) override;
  Real computeQpOffDiagJacobianScalar(unsigned int jvar) override;
  const unsigned int _component;
  const unsigned int _displacement_number;
  const unsigned int _rotation_number;
  const VariableValue & _displacement;
  const VariableValue & _rotations;
  const libMesh::Point _reference_point;
  const Function & _translation;
};
