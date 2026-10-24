#pragma once
#include "NodalConstraint.h"
#include "MeshGenerator.h"
#include "AqueductMPCAssemblyStorage.h"
#include <memory>

// Mixed linear row C*u-compliance*lambda=0 and conjugate C^T*lambda; zero compliance is kinematic.
class AqueductLinearMPC : public NodalConstraint
{
public:
  static InputParameters validParams();
  AqueductLinearMPC(const InputParameters &);
  void computeResidual(const NumericVector<Number> &) override;
  void computeJacobian(const SparseMatrix<Number> &) override;
protected:
  Real computeQpResidual(Moose::ConstraintType) override;
  Real computeQpJacobian(Moose::ConstraintJacobianType) override;
  std::vector<dof_id_type> dofs();
  void dofs(std::vector<dof_id_type> &);
  const dof_id_type _anchor;
  std::vector<dof_id_type> _nodes;
  std::vector<Real> _coefficients;
  std::vector<MooseVariable *> _variables;
  MooseVariable * _multiplier;
  const Real _compliance;
  std::unique_ptr<AqueductMPCAssemblyStorage> _storage;
  std::vector<Real> _scales_scratch;
};

// Zero-dimensional support activates multiplier DOFs only at actual anchors.
class AqueductConstraintNodesGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();
  AqueductConstraintNodesGenerator(const InputParameters &);
  std::unique_ptr<MeshBase> generate() override;
protected:
  std::unique_ptr<MeshBase> & _input;
};
