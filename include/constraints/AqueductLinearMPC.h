#pragma once
#include "NodalConstraint.h"
#include "MeshGenerator.h"

// Exact linear constraint C*u=0 with its conjugate nodal multiplier: no penalty.
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
  const dof_id_type _anchor;
  std::vector<dof_id_type> _nodes;
  std::vector<Real> _coefficients;
  std::vector<MooseVariable *> _variables;
  MooseVariable * _multiplier;
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
