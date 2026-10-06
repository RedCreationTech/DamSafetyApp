#include "RigidPlaneMomentConstraint.h"
#include "DelimitedFileReader.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include "MooseVariableScalar.h"
#include "SystemBase.h"
#include "libmesh/numeric_vector.h"
#include <map>
#include <cmath>

registerMooseObject("DamSafetyApp", RigidPlaneMomentConstraint);

InputParameters RigidPlaneMomentConstraint::validParams()
{
  auto params = NodalScalarKernel::validParams();
  params.addClassDescription("Moment equilibrium for three free rotations of a coupled rigid plane.");
  params.addRequiredCoupledVar("tractions", "Three nodal LM traction fields in x/y/z order.");
  params.addRequiredParam<FileName>("moment_file", "Reference surface nodal-integral CSV.");
  params.addParam<SubdomainID>("surface_block", 500, "Lower-dimensional RP surface block.");
  // This replicated-mesh wall also has non-topological embedded node/element
  // hosts. Keep their solution DOFs available on every rank before upstream
  // embedded constraints reinitialize neighbor variables. Algebraic ghosting
  // does not create a dense matrix coupling graph or change the equations.
  params.addRelationshipManager("GhostEverything", Moose::RelationshipManagerType::ALGEBRAIC);
  return params;
}

RigidPlaneMomentConstraint::RigidPlaneMomentConstraint(const InputParameters & parameters)
  : NodalScalarKernel(parameters), _traction_numbers(coupledIndices("tractions")),
    _surface_block(getParam<SubdomainID>("surface_block"))
{
  if (_traction_numbers.size() != 3)
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
  // Source labels are evidence keys, not persistent runtime mesh node IDs.
  for (const auto & entry : data)
  {
    _reference_nodes.emplace_back(entry.second[0], entry.second[1], entry.second[2]);
    _moments.emplace_back(entry.second[3], entry.second[4], entry.second[5]);
  }
}

void RigidPlaneMomentConstraint::initialSetup()
{
  std::set<const Node *> nodes;
  for (const auto * elem : _mesh.getMesh().active_element_ptr_range())
    if (elem->subdomain_id() == _surface_block)
      for (unsigned int n = 0; n < elem->n_nodes(); ++n)
        nodes.insert(elem->node_ptr(n));
  if (nodes.size() != _reference_nodes.size())
    mooseError("RP surface has ", nodes.size(), " nodes, expected ", _reference_nodes.size());
  for (const auto & point : _reference_nodes)
  {
    const Node * match = nullptr;
    for (const auto * node : nodes)
      if ((*node - point).norm() <= 1e-10)
      {
        if (match) mooseError("Nonunique RP surface coordinate mapping.");
        match = node;
      }
    if (!match) mooseError("Original RP surface point is absent: ", point);
    _surface_nodes.push_back(match);
  }
  for (unsigned int d = 0; d < 3; ++d)
    for (unsigned int n = 0; n < _surface_nodes.size(); ++n)
      tractionDof(d, n);
  mooseInfo("RP surface coordinate identity verified for ", _surface_nodes.size(), " nodes.");
}

dof_id_type RigidPlaneMomentConstraint::tractionDof(unsigned int traction, unsigned int n) const
{
  const auto * node = _surface_nodes.at(n);
  if (!node || node->n_dofs(_sys.number(), _traction_numbers.at(traction)) != 1)
    mooseError("RP surface node ", node->id(), " at ", *node,
               " must have exactly one traction DOF in system ", _sys.number(),
               " variable ", _traction_numbers.at(traction));
  return node->dof_number(_sys.number(), _traction_numbers[traction], 0);
}

void RigidPlaneMomentConstraint::reinit()
{
  // Scalar-kernel reinit is called on every rank before the scalar-DOF owner
  // assembles its rows. The surface moments need all 224 nodal tractions;
  // geometric mesh replication does not replicate a distributed solution.
  // Collect collectively here, never inside owner-only computeResidual().
  if (_communicator.size() > 1)
    _sys.currentSolution()->localize(_parallel_solution);
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
  if (_local_re.size() != 3)
    mooseError("RP moment residual must have three rows, got ", _local_re.size());
  for (unsigned int r = 0; r < 3; ++r)
  {
    _local_re(r) = 0.0;
    for (unsigned int d = 0; d < 3; ++d)
      for (unsigned int n = 0; n < _node_ids.size(); ++n)
      {
        const auto dof = tractionDof(d, n);
        const auto traction = _communicator.size() > 1 ? _parallel_solution.at(dof)
                                                       : (*_sys.currentSolution())(dof);
        _local_re(r) += coefficient(r, d, n) * traction;
      }
  }
  assignTaggedLocalResidual();
}

void RigidPlaneMomentConstraint::computeJacobian()
{
  const auto & rows = _var.dofIndices();
  if (rows.size() != 3)
    mooseError("RP moment Jacobian must have three scalar DOFs.");
  // These are global surface moments, not an element-local field block.
  // Cache entries using the actual nodal global DOFs, avoiding stale/local
  // block dimensions and the finite capacity of a single global AD vector.
  for (unsigned int d = 0; d < 3; ++d)
  {
    for (unsigned int r = 0; r < 3; ++r)
      for (unsigned int n = 0; n < _node_ids.size(); ++n)
      {
        const Real value = coefficient(r, d, n);
        if (value != 0)
          addJacobianElement(_assembly, value, rows[r], tractionDof(d, n), _var.scalingFactor());
      }
  }
}
