#include "ComputeUniformHex8IncrementalStrain.h"
#include "MooseVariable.h"
#include "libmesh/quadrature.h"

registerMooseObject("DamSafetyApp", ComputeUniformHex8IncrementalStrain);

InputParameters ComputeUniformHex8IncrementalStrain::validParams()
{
  auto params = ComputeIncrementalStrain::validParams();
  params.addClassDescription("Reference-geometry uniform-strain HEX8, one material point.");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}

ComputeUniformHex8IncrementalStrain::ComputeUniformHex8IncrementalStrain(
    const InputParameters & parameters) : ComputeIncrementalStrain(parameters)
{
  if (_ndisp != 3 || _volumetric_locking_correction || getParam<bool>("use_displaced_mesh"))
    mooseError("Uniform HEX8 requires three displacements, reference geometry and no B-bar.");
  for (unsigned int d = 0; d < 3; ++d)
    _nodal_displacements[d] = getVar("displacements", d);
}

void ComputeUniformHex8IncrementalStrain::computeTotalStrainIncrement(RankTwoTensor & increment)
{
  if (_current_elem->type() != libMesh::HEX8 || _qrule->n_points() != 1)
    mooseError("Uniform HEX8 requires HEX8 topology and exactly one constitutive point.");
  auto it = _operators.find(_current_elem->id());
  if (it == _operators.end())
  {
    std::array<libMesh::Point, 8> points;
    for (unsigned int i = 0; i < 8; ++i)
      points[i] = _current_elem->point(i);
    it = _operators.emplace(_current_elem->id(), UniformHex8::compute(points)).first;
  }
  RankTwoTensor gradient, old_gradient;
  for (unsigned int d = 0; d < 3; ++d)
    for (unsigned int i = 0; i < 8; ++i)
    {
      const auto & node = _current_elem->node_ref(i);
      const Real current = _nodal_displacements[d]->getNodalValue(node);
      const Real old = _nodal_displacements[d]->getNodalValueOld(node);
      for (unsigned int a = 0; a < 3; ++a)
      {
        gradient(d, a) += current * it->second.gradient[i](a);
        old_gradient(d, a) += old * it->second.gradient[i](a);
      }
    }
  _deformation_gradient[_qp] = gradient;
  _deformation_gradient[_qp].addIa(1.0);
  const auto delta = gradient - old_gradient;
  increment = (delta + delta.transpose()) * 0.5;
}
