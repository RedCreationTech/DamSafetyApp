#pragma once
#include "NodalConstraint.h"
#include "AqueductMPCAssemblyStorage.h"
#include "AqueductMPCCrossJacobian.h"
#include <array>

/** Audited host interpolation with explicit conjugate multipliers.
 * Physical mass/force remain in their original tags. All algebraic rows are
 * current-time equations, so the whole-system HHT can preserve physical forces.
 */
class AqueductEmbeddedMPC : public NodalConstraint
{
public:
  static InputParameters validParams();
  AqueductEmbeddedMPC(const InputParameters &);
  void computeResidual(const NumericVector<Number> &) override;
  void computeJacobian(const SparseMatrix<Number> &) override;
protected:
  struct Row
  {
    std::array<dof_id_type,9> nodes;
    std::array<Real,9> coefficients;
  };
  std::vector<Row> _rows;
  const bool _reuse_storage;
  std::vector<AqueductMPCAssemblyStorage> _storage;
  std::vector<AqueductMPCCrossJacobian> _cross_jacobians;
  std::array<MooseVariable *,3> _physical, _multipliers;
  std::vector<dof_id_type> dofs(const Row &,unsigned int) const;
  void dofs(const Row &,unsigned int,std::vector<dof_id_type> &) const;
  Real computeQpResidual(Moose::ConstraintType) override;
  Real computeQpJacobian(Moose::ConstraintJacobianType) override;
};
