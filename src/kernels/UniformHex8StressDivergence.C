#include "UniformHex8StressDivergence.h"
#include "MooseVariable.h"
#include "libmesh/quadrature.h"

registerMooseObject("DamSafetyApp", UniformHex8StressDivergence);

InputParameters UniformHex8StressDivergence::validParams()
{
  auto params = StressDivergenceTensors::validParams();
  params.addClassDescription("Small-strain uniform HEX8 forces and total-stiffness hourglass control.");
  params.addRequiredRangeCheckedParam<Real>("initial_shear_modulus", "initial_shear_modulus > 0",
                                           "Initial elastic shear modulus in Pa, retained under damage.");
  params.addRangeCheckedParam<Real>("hourglass_factor", 0.005, "hourglass_factor >= 0",
                                   "Total-stiffness hourglass factor; 0.005 is the Standard default.");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

UniformHex8StressDivergence::UniformHex8StressDivergence(const InputParameters & parameters)
  : StressDivergenceTensors(parameters),
    _initial_shear_modulus(getParam<Real>("initial_shear_modulus")),
    _hourglass_factor(getParam<Real>("hourglass_factor")),
    _nodal_displacement(getVar("displacements", _component))
{
  if (_ndisp != 3 || _volumetric_locking_correction || _use_displaced_mesh ||
      _use_finite_deform_jacobian)
    mooseError("Uniform HEX8 supports only three-component reference small strain, without B-bar.");
}

const UniformHex8::Operators & UniformHex8StressDivergence::operators()
{
  if (_current_elem->type() != libMesh::HEX8 || _qrule->n_points() != 1 ||
      _test.size() != 8 || _phi.size() != 8)
    mooseError("Uniform HEX8 requires first-order HEX8 and one constitutive integration point.");
  auto it = _operators.find(_current_elem->id());
  if (it == _operators.end())
  {
    std::array<libMesh::Point, 8> points;
    for (unsigned int i = 0; i < 8; ++i)
      points[i] = _current_elem->point(i);
    it = _operators.emplace(_current_elem->id(), UniformHex8::compute(points)).first;
  }
  return it->second;
}

Real UniformHex8StressDivergence::computeQpResidual()
{
  const auto & op = operators();
  Real force = 0.0;
  for (unsigned int d = 0; d < 3; ++d)
    force += op.volume * _stress[_qp](_component, d) * op.gradient[_i](d);
  for (unsigned int j = 0; j < 8; ++j)
    force += op.hourglassStiffness(_i, j, _initial_shear_modulus, _hourglass_factor) *
             _nodal_displacement->getNodalValue(_current_elem->node_ref(j));
  // Kernel assembly multiplies by the centroid quadrature weight. Undo that
  // weight: uniform strain uses the true integrated reference volume above.
  return force / (_JxW[_qp] * _coord[_qp]);
}

Real UniformHex8StressDivergence::stiffness(unsigned int column)
{
  const auto & op = operators();
  Real value = 0.0;
  for (unsigned int k = 0; k < 3; ++k)
    for (unsigned int l = 0; l < 3; ++l)
      value += op.volume * op.gradient[_i](k) *
               _Jacobian_mult[_qp](_component, k, column, l) * op.gradient[_j](l);
  if (column == _component)
    value += op.hourglassStiffness(_i, _j, _initial_shear_modulus, _hourglass_factor);
  return value / (_JxW[_qp] * _coord[_qp]);
}

Real UniformHex8StressDivergence::computeQpJacobian() { return stiffness(_component); }

Real UniformHex8StressDivergence::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int d = 0; d < 3; ++d)
    if (jvar == _disp_var[d])
      return stiffness(d);
  return 0.0;
}
