#pragma once

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

/** Constant mixed-row storage. DOF values must still be refreshed on every assembly. */
class AqueductMPCAssemblyStorage
{
public:
  explicit AqueductMPCAssemblyStorage(std::vector<libMesh::Real> coefficients,
                                      bool assemble_compliance = true)
    : _coefficients(std::move(coefficients)),
      _assemble_compliance(assemble_compliance),
      _scales(_coefficients.size(), std::numeric_limits<libMesh::Real>::quiet_NaN()),
      _residual(_coefficients.size() + 1),
      _jacobian(_coefficients.size() + 1, _coefficients.size() + 1)
  {
    if (_coefficients.empty())
      throw std::invalid_argument("MPC storage requires a nonempty row");
    _dofs.reserve(_coefficients.size() + 1);
  }

  std::vector<libMesh::dof_id_type> & dofs()
  {
    _dofs.clear();
    return _dofs;
  }

  libMesh::DenseVector<libMesh::Number> & residual()
  {
    _residual.zero();
    return _residual;
  }

  template <typename Scales>
  const libMesh::DenseMatrix<libMesh::Number> &
  jacobian(const Scales & scales, libMesh::Real multiplier_scale, libMesh::Real compliance)
  {
    if (scales.size() != _coefficients.size())
      throw std::invalid_argument("MPC storage scaling size differs from row");
    bool changed = !_valid || multiplier_scale != _multiplier_scale || compliance != _compliance;
    for (std::size_t j = 0; j < scales.size(); ++j)
      changed = changed || scales[j] != _scales[j];
    if (changed)
    {
      const auto n = _coefficients.size();
      _jacobian.zero();
      for (std::size_t j = 0; j < n; ++j)
      {
        _jacobian(j, n) = _coefficients[j] * scales[j];
        _jacobian(n, j) = _coefficients[j] * multiplier_scale;
        _scales[j] = scales[j];
      }
      if (_assemble_compliance)
        _jacobian(n, n) = -compliance * multiplier_scale;
      _multiplier_scale = multiplier_scale;
      _compliance = compliance;
      _valid = true;
      ++_rebuilds;
    }
    return _jacobian;
  }

  std::size_t rebuilds() const { return _rebuilds; }

private:
  const std::vector<libMesh::Real> _coefficients;
  const bool _assemble_compliance;
  std::vector<libMesh::Real> _scales;
  std::vector<libMesh::dof_id_type> _dofs;
  libMesh::DenseVector<libMesh::Number> _residual;
  libMesh::DenseMatrix<libMesh::Number> _jacobian;
  libMesh::Real _multiplier_scale = 0., _compliance = 0.;
  bool _valid = false;
  std::size_t _rebuilds = 0;
};
