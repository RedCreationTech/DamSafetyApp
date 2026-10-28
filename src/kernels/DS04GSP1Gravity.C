#include "DS04GSP1Gravity.h"
#include "ActionWarehouse.h"
#include "FinitePlaneStrainNewmarkPhysics.h"
#include "Function.h"
#include "MooseApp.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "libmesh/elem.h"

registerMooseObject("DamSafetyApp", DS04GSP1Gravity);

InputParameters DS04GSP1Gravity::validParams()
{
  auto p = Gravity::validParams();
  p.addClassDescription("SP1-only optional current-area gravity Jacobian; original Gravity residual");
  p.addRequiredCoupledVar("displacements", "The same two FIRST LAGRANGE SP1 displacements");
  return p;
}

DS04GSP1Gravity::DS04GSP1Gravity(const InputParameters & p)
  : Gravity(p), _displacement_vars{{coupled("displacements", 0), coupled("displacements", 1)}}
{
  if (coupledComponents("displacements") != 2 || _mesh.dimension() != 2 ||
      !_use_displaced_mesh || _var.number() != _displacement_vars[1])
    mooseError("SP1 gravity requires current-geometry Y force and two XY displacements");
  for (unsigned b = 0; b < 2; ++b)
    if (getVar("displacements", b)->feType() != _var.feType() ||
        _var.feType().family != libMesh::LAGRANGE || _var.feType().order != libMesh::FIRST)
      mooseError("SP1 gravity requires matching FIRST LAGRANGE displacement spaces");
}

void DS04GSP1Gravity::initialSetup()
{
  Gravity::initialSetup();
  unsigned count = 0;
  for (const auto * physics : _app.actionWarehouse().getActions<FinitePlaneStrainNewmarkPhysics>())
  {
    ++count;
    if (physics->getParam<MooseEnum>("spatial_discretization") != "sp1_center_hg")
      mooseError("DS04GSP1Gravity requires the SP1 physics; LEGACY is unchanged");
    _complete = physics->getParam<bool>("sp1_gravity_jacobian");
  }
  if (count != 1 || getBlockCoordSystem() != Moose::COORD_XYZ)
    mooseError("SP1 gravity requires one Cartesian SP1 physics");
}

Real DS04GSP1Gravity::tangent(unsigned direction)
{
  if (!_complete)
    return 0;
  if (_current_elem->type() != QUAD4 || _qrule->n_points() != 1)
    mooseError("SP1 gravity Jacobian requires one center QP on QUAD4");
  // With uniform constant gravity/density, only the current volume measure varies.
  return _density[_qp] * _test[_i][_qp] * -_value *
         _function.value(_t + _alpha * _dt, _q_point[_qp]) * _grad_phi[_j][_qp](direction);
}

Real DS04GSP1Gravity::computeQpJacobian() { return tangent(1); }

Real DS04GSP1Gravity::computeQpOffDiagJacobian(unsigned jvar)
{
  for (unsigned b = 0; b < 2; ++b)
    if (jvar == _displacement_vars[b])
      return tangent(b);
  return 0;
}
