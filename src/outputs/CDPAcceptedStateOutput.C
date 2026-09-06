#include "CDPAcceptedStateOutput.h"
#include "CDPAssemblyProbe.h"
#include "FEProblemBase.h"
#include "MaterialPropertyStorage.h"
#include "MooseMesh.h"
#include "NonlinearSystemBase.h"
#include "RankTwoTensor.h"
#include "MooseVariableFieldBase.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include <algorithm>
#include <fstream>
#include <iomanip>

registerMooseObject("DamSafetyApp", CDPAcceptedStateOutput);
InputParameters CDPAcceptedStateOutput::validParams()
{
  auto p = FileOutput::validParams();
  p.set<ExecFlagEnum>("execute_on", true) = EXEC_TIMESTEP_END;
  p.addParam<std::vector<Real>>("fd_targets", {0., 0., 0.11, 2., 0., 0., 0.14, 2.},
      "Flat x y z component tuples for diagnostic columns; component 0/1/2 is X/Y/Z.");
  p.addParam<std::vector<Real>>("capture_times", {}, "Optional exact accepted times; never forces steps");
  p.addParam<bool>("assembly_probe", false, "Capture passive assembly material data");
  p.addClassDescription("Capture accepted global solution, assembled residual/Jacobian and complete stored material histories.");
  return p;
}
CDPAcceptedStateOutput::CDPAcceptedStateOutput(const InputParameters & p) : FileOutput(p) {}
std::string CDPAcceptedStateOutput::filename()
{
  return _file_base + "_step" + std::to_string(_t_step);
}
void CDPAcceptedStateOutput::writeHistory(const std::string & stage)
{
  const auto & storage = _problem_ptr->getMaterialPropertyStorage();
  std::ofstream f(filename() + "_" + stage + "_rank" + std::to_string(processor_id()) + ".csv");
  f << std::setprecision(17) << "time,dt,step,element,side,state,property,qp,component,value\n";
  for (auto state : storage.stateIndexRange())
    for (const auto & entry : storage.props(state))
      for (const auto & side : entry.second)
        for (unsigned int id = 0; id < side.second.size(); ++id)
        {
          const auto * value = side.second.queryValue(id);
          if (!value) continue;
          const auto name = storage.queryStatefulPropName(storage.statefulProps().at(id));
          if (!name) mooseError("Missing stateful property name");
          const auto * scalar = dynamic_cast<const MaterialProperty<Real> *>(value);
          const auto * tensor = dynamic_cast<const MaterialProperty<RankTwoTensor> *>(value);
          if (!scalar && !tensor) mooseError("Unsupported diagnostic state type: ", *name);
          for (unsigned int qp=0; qp<value->size(); ++qp)
            for (unsigned int c=0; c<(scalar ? 1u : 9u); ++c)
              f << _time << ',' << _dt << ',' << _t_step << ',' << entry.first->id() << ','
                << side.first << ',' << state << ',' << *name << ',' << qp << ',' << c << ','
                << (scalar ? (*scalar)[qp] : (*tensor)[qp](c/3,c%3)) << '\n';
        }
  if (!f) mooseError("Cannot write accepted material state");
}
void CDPAcceptedStateOutput::output()
{
  const auto & times = getParam<std::vector<Real>>("capture_times");
  if (!times.empty() && std::none_of(times.begin(), times.end(),
      [this](Real t) { return std::abs(_time-t) < 1e-9; })) return;
  const bool probe = getParam<bool>("assembly_probe");
  const std::string rank = "_rank" + std::to_string(processor_id()) + ".csv";
  auto & nl = _problem_ptr->getNonlinearSystemBase(0);
  auto & sys = dynamic_cast<libMesh::NonlinearImplicitSystem &>(nl.system());
  const auto stem = filename();
  // Element assembly needs the ghosted current solution, not the distributed
  // owning System::solution vector (remote element DOFs are otherwise unavailable).
  const auto * accepted = nl.currentSolution();
  auto solution = accepted->clone();
  solution->print_matlab(stem + "_u.m");
  nl.solutionOld().print_matlab(stem + "_u_old.m");
  nl.solutionOlder().print_matlab(stem + "_u_older.m");
  std::ofstream map(stem + "_map_rank" + std::to_string(processor_id()) + ".csv");
  map << std::setprecision(17) << "node,x,y,z,variable,dof,scaling\n";
  for (const auto * node : _mesh_ptr->getMesh().local_node_ptr_range())
    for (unsigned int v=0; v<sys.n_vars(); ++v)
      if (node->n_dofs(sys.number(),v))
        map << node->id() << ',' << (*node)(0) << ',' << (*node)(1) << ',' << (*node)(2)
            << ',' << sys.variable_name(v) << ',' << node->dof_number(sys.number(),v,0) << ',' << nl.getVariable(0,v).scalingFactor() << '\n';
  std::ofstream qpfile(stem + "_qp_rank" + std::to_string(processor_id()) + ".csv");
  qpfile << std::setprecision(17) << "element,qp,x,y,z\n";
  libMesh::QGauss q(3,libMesh::SECOND);
  auto fe = libMesh::FEBase::build(3,libMesh::FEType(libMesh::FIRST,libMesh::LAGRANGE));
  fe->attach_quadrature_rule(&q);
  const auto & xyz = fe->get_xyz();
  for (const auto * elem : _mesh_ptr->getMesh().active_local_element_ptr_range())
  {
    fe->reinit(elem);
    for (unsigned int i=0; i<xyz.size(); ++i)
      qpfile << elem->id() << ',' << i << ',' << xyz[i](0) << ',' << xyz[i](1) << ',' << xyz[i](2) << '\n';
  }
  writeHistory("before");
  auto residual = solution->zero_clone();
  auto jacobian = sys.get_system_matrix().zero_clone();
  _problem_ptr->computeResidual(*accepted,*residual,0);
  residual->print_matlab(stem + "_residual.m");
  if (probe) CDPAssemblyProbe::beginCapture(stem + "_ip_base" + rank, true);
  _problem_ptr->computeJacobian(*accepted,*jacobian,0);
  if (probe) CDPAssemblyProbe::endCapture();
  jacobian->print_matlab(stem + "_jacobian.m");
  // Configurable physical coordinates avoid reliance on MPI-dependent DOF IDs.
  // Keep the original two centerline Z columns as the default.
  const auto & targets = getParam<std::vector<Real>>("fd_targets");
  if (targets.empty() || targets.size() % 4)
    paramError("fd_targets", "Expected nonempty x y z component tuples");
  for (unsigned int target = 0; target < targets.size(); target += 4)
  {
    const Real component = targets[target + 3];
    if (component != 0. && component != 1. && component != 2.)
      paramError("fd_targets", "Component must be 0, 1 or 2");
    const std::string variable = component == 0. ? "disp_x" : component == 1. ? "disp_y" : "disp_z";
    libMesh::dof_id_type dof = libMesh::DofObject::invalid_id;
    for (const auto * node : _mesh_ptr->getMesh().local_node_ptr_range())
      if (std::abs((*node)(0)-targets[target]) < 1e-12 && std::abs((*node)(1)-targets[target+1]) < 1e-12 && std::abs((*node)(2)-targets[target+2]) < 1e-12)
        dof = node->dof_number(sys.number(),sys.variable_number(variable),0);
    _communicator.min(dof);
    if (dof == libMesh::DofObject::invalid_id) mooseError("Missing mirror diagnostic node");
    for (const Real h : {1e-10, 1e-11})
    {
      auto plus = solution->clone();
      auto minus = solution->clone();
      if (dof >= plus->first_local_index() && dof < plus->last_local_index())
      {
        plus->add(dof,h);
        minus->add(dof,-h);
      }
      plus->close(); minus->close();
      auto rp = residual->zero_clone();
      auto rm = residual->zero_clone();
      const std::string ip_stem = stem + "_ip_dof" + std::to_string(dof) + (h==1e-10 ? "_h10" : "_h11");
      if (probe) CDPAssemblyProbe::beginCapture(ip_stem + "_plus" + rank, false);
      _problem_ptr->computeResidual(*plus,*rp,0);
      if (probe) CDPAssemblyProbe::endCapture();
      if (probe) CDPAssemblyProbe::beginCapture(ip_stem + "_minus" + rank, false);
      _problem_ptr->computeResidual(*minus,*rm,0);
      if (probe) CDPAssemblyProbe::endCapture();
      rp->add(-1.,*rm); rp->scale(0.5/h);
      rp->print_matlab(stem + "_fd_dof" + std::to_string(dof) + (h==1e-10 ? "_h10.m" : "_h11.m"));
      // Reset the problem's current solution pointer before scratch vectors die.
      _problem_ptr->computeResidual(*accepted,*residual,0);
    }
  }
  // Re-evaluate the original accepted solution last, with original old histories.
  _problem_ptr->computeResidual(*accepted,*residual,0);
  writeHistory("after");
}
