#include "RigidPlaneConstraint.h"

registerMooseObject("DamSafetyApp", RigidPlaneConstraint);

InputParameters RigidPlaneConstraint::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription("Small-rotation rigid-plane kinematic equation with free RP rotations.");
  params.addRequiredRangeCheckedParam<unsigned int>("component", "component < 3", "Displacement axis.");
  params.addRequiredCoupledVar("displacement", "The displacement component being constrained.");
  params.addRequiredCoupledVar("rotations", "A THIRD-order SCALAR variable with three RP rotations.");
  params.addRequiredParam<libMesh::Point>("reference_point", "Original RP position in metres.");
  params.addParam<FunctionName>("translation", "0", "Prescribed RP displacement component in metres.");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

RigidPlaneConstraint::RigidPlaneConstraint(const InputParameters & parameters)
  : Kernel(parameters), _component(getParam<unsigned int>("component")),
    _displacement_number(coupled("displacement")), _rotation_number(coupledScalar("rotations")),
    _displacement(coupledValue("displacement")), _rotations(coupledScalarValue("rotations")),
    _reference_point(getParam<libMesh::Point>("reference_point")),
    _translation(getFunction("translation"))
{
  if (getParam<bool>("use_displaced_mesh"))
    paramError("use_displaced_mesh", "Rigid plane uses reference small-strain geometry.");
}

Real RigidPlaneConstraint::computeQpResidual()
{
  if (_rotations.size() != 3 || _current_elem->dim() != 2)
    mooseError("RigidPlaneConstraint requires three scalar rotations and a surface block.");
  const auto r = _q_point[_qp] - _reference_point;
  const unsigned int a = (_component + 1) % 3;
  const unsigned int b = (_component + 2) % 3;
  const Real rotation_displacement = _rotations[a] * r(b) - _rotations[b] * r(a);
  return _test[_i][_qp] * (_displacement[_qp] - _translation.value(_t, _q_point[_qp]) -
                         rotation_displacement);
}

Real RigidPlaneConstraint::computeQpOffDiagJacobian(unsigned int jvar)
{
  return jvar == _displacement_number ? _test[_i][_qp] * _phi[_j][_qp] : 0.0;
}

Real RigidPlaneConstraint::computeQpOffDiagJacobianScalar(unsigned int jvar)
{
  if (jvar != _rotation_number) return 0.0;
  const auto r = _q_point[_qp] - _reference_point;
  const unsigned int a = (_component + 1) % 3;
  const unsigned int b = (_component + 2) % 3;
  const Real coefficient = _j == a ? -r(b) : _j == b ? r(a) : 0.0;
  return _test[_i][_qp] * coefficient;
}
