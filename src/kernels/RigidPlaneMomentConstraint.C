#include "RigidPlaneMomentConstraint.h"
#include "DelimitedFileReader.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include <map>
#include <cmath>

registerMooseObject("DamSafetyApp", RigidPlaneMomentConstraint);

InputParameters RigidPlaneMomentConstraint::validParams()
{
  auto params = NodalScalarKernel::validParams();
  params.addClassDescription("Moment equilibrium for three free rotations of a coupled rigid plane.");
  params.addRequiredCoupledVar("tractions", "Three nodal LM traction fields in x/y/z order.");
  params.addRequiredParam<FileName>("moment_file", "Reference surface nodal-integral CSV.");
  return params;
}

RigidPlaneMomentConstraint::RigidPlaneMomentConstraint(const InputParameters & parameters)
  : NodalScalarKernel(parameters), _traction_numbers(coupledIndices("tractions")),
    _tractions(coupledValues("tractions"))
{
  if (_tractions.size() != 3)
    paramError("tractions", "Supply three traction fields in x/y/z order.");
  MooseUtils::DelimitedFileReader reader(getParam<FileName>("moment_file"));
  reader.setHeaderFlag(MooseUtils::DelimitedFileReader::HeaderFlag::ON);
  reader.read();
  const auto & ids = reader.getData("moose_node_id");
  const auto & x = reader.getData("x"); const auto & y = reader.getData("y");
  const auto & z = reader.getData("z");
  const auto & mx = reader.getData("moment_x"); const auto & my = reader.getData("moment_y");
  const auto & mz = reader.getData("moment_z");
  std::map<dof_id_type, std::array<Real, 6>> data;
  for (unsigned int i = 0; i < ids.size(); ++i)
  {
    if (ids[i] < 0 || ids[i] != std::floor(ids[i]))
      paramError("moment_file", "Invalid node ID.");
    if (!data.emplace(ids[i], std::array<Real,6>{x[i],y[i],z[i],mx[i],my[i],mz[i]}).second)
      paramError("moment_file", "Duplicate node ID.");
  }
  if (data.size() != _node_ids.size())
    paramError("moment_file", "Moment rows must match the complete top nodeset.");
  for (const auto id : _node_ids)
  {
    const auto it = data.find(id);
    const auto * node = _mesh.getMesh().query_node_ptr(id);
    if (it == data.end() || !node)
      paramError("moment_file", "Unknown top node ", id);
    for (unsigned int d = 0; d < 3; ++d)
      if (std::abs((*node)(d) - it->second[d]) > 1e-10)
        paramError("moment_file", "Node numbering/coordinates do not match the input mesh.");
    _moments.emplace_back(it->second[3], it->second[4], it->second[5]);
  }
}

Real RigidPlaneMomentConstraint::coefficient(unsigned int rotation, unsigned int traction,
                                           unsigned int node) const
{
  const unsigned int a = (rotation + 1) % 3;
  const unsigned int b = (rotation + 2) % 3;
  return traction == a ? _moments[node](b) : traction == b ? -_moments[node](a) : 0.0;
}

void RigidPlaneMomentConstraint::computeResidual()
{
  if (_u.size() != 3) mooseError("RP rotations must be a THIRD-order SCALAR variable.");
  prepareVectorTag(_assembly, _var.number());
  for (unsigned int r = 0; r < 3; ++r)
  {
    _local_re(r) = 0.0;
    for (unsigned int d = 0; d < 3; ++d)
      for (unsigned int n = 0; n < _node_ids.size(); ++n)
        _local_re(r) += coefficient(r, d, n) * (*_tractions[d])[n];
  }
  assignTaggedLocalResidual();
}

void RigidPlaneMomentConstraint::computeJacobian()
{
  prepareMatrixTag(_assembly, _var.number(), _var.number());
  _local_ke.zero();
  assignTaggedLocalMatrix();
  for (unsigned int d = 0; d < 3; ++d)
  {
    prepareMatrixTag(_assembly, _var.number(), _traction_numbers[d]);
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int n = 0; n < _node_ids.size(); ++n)
        _local_ke(r,n) = coefficient(r,d,n);
    assignTaggedLocalMatrix();
  }
}
