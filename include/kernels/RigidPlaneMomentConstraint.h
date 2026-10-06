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
  void initialSetup() override;
  void reinit() override {}
  void computeResidual() override;
  void computeJacobian() override;
protected:
  Real coefficient(unsigned int rotation, unsigned int traction, unsigned int node) const;
  dof_id_type tractionDof(unsigned int traction, unsigned int node) const;
  const std::vector<unsigned int> _traction_numbers;
  const SubdomainID _surface_block;
  std::vector<Point> _reference_nodes;
  std::vector<const Node *> _surface_nodes;
  std::vector<RealVectorValue> _moments;
};
