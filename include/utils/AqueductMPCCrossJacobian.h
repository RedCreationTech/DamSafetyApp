#pragma once

#include "libmesh/dense_matrix.h"
#include "libmesh/dof_map.h"
#include <algorithm>
#include <stdexcept>
#include <vector>

/** Fixed mixed row only: the physical/physical block is identically zero. */
class AqueductMPCCrossJacobian
{
public:
  explicit AqueductMPCCrossJacobian(std::size_t physical_terms)
    : _n(physical_terms),
      _physical_multiplier(_n, 1),
      _multiplier_physical(1, _n),
      _compliance(1, 1),
      _multiplier_dof(1)
  {
    if (!_n)
      throw std::invalid_argument("MPC cross blocks require physical terms");
    _physical_dofs.reserve(_n);
  }

  template <typename Assemble>
  void assemble(const libMesh::DenseMatrix<libMesh::Number> & matrix,
                const std::vector<libMesh::dof_id_type> & dofs,
                const libMesh::DofMap & dof_map,
                Assemble add_block)
  {
    if (matrix.m() != _n + 1 || matrix.n() != _n + 1 || dofs.size() != _n + 1)
      throw std::invalid_argument("MPC cross block dimensions differ from the fixed row");

    // General libMesh constraints can create extra entries: preserve the complete old path.
    if (std::any_of(dofs.begin(), dofs.end(),
                    [&](auto id) { return dof_map.is_constrained_dof(id); }))
    {
      add_block(matrix, dofs, dofs);
      return;
    }

    // Refresh values after every DOF reinitialization; keep the old nonzero insertion order.
    _physical_dofs.assign(dofs.begin(), dofs.end() - 1);
    _multiplier_dof[0] = dofs.back();
    for (std::size_t j = 0; j < _n; ++j)
    {
      _physical_multiplier(j, 0) = matrix(j, _n);
      _multiplier_physical(0, j) = matrix(_n, j);
    }
    add_block(_physical_multiplier, _physical_dofs, _multiplier_dof);
    add_block(_multiplier_physical, _multiplier_dof, _physical_dofs);
    if (matrix(_n, _n) != 0.)
    {
      _compliance(0, 0) = matrix(_n, _n);
      add_block(_compliance, _multiplier_dof, _multiplier_dof);
    }
  }

private:
  const std::size_t _n;
  libMesh::DenseMatrix<libMesh::Number> _physical_multiplier, _multiplier_physical, _compliance;
  std::vector<libMesh::dof_id_type> _physical_dofs, _multiplier_dof;
};
