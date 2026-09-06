#include "CDPAcceptedStateOutput.h"
#include "FEProblemBase.h"
#include "MaterialPropertyStorage.h"
#include "MooseMesh.h"
#include "NonlinearSystemBase.h"
#include "RankTwoTensor.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include <fstream>
#include <iomanip>

registerMooseObject("DamSafetyApp", CDPAcceptedStateOutput);
InputParameters CDPAcceptedStateOutput::validParams()
{
  auto p = FileOutput::validParams();
  p.set<ExecFlagEnum>("execute_on", true) = EXEC_TIMESTEP_END;
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
  auto & nl = _problem_ptr->getNonlinearSystemBase(0);
  auto & sys = dynamic_cast<libMesh::NonlinearImplicitSystem &>(nl.system());
  const auto stem = filename();
  auto solution = sys.solution->clone();
  solution->print_matlab(stem + "_u.m");
  nl.solutionOld().print_matlab(stem + "_u_old.m");
  nl.solutionOlder().print_matlab(stem + "_u_older.m");
  std::ofstream map(stem + "_map_rank" + std::to_string(processor_id()) + ".csv");
  map << std::setprecision(17) << "node,x,y,z,variable,dof\n";
  for (const auto * node : _mesh_ptr->getMesh().local_node_ptr_range())
    for (unsigned int v=0; v<sys.n_vars(); ++v)
      if (node->n_dofs(sys.number(),v))
        map << node->id() << ',' << (*node)(0) << ',' << (*node)(1) << ',' << (*node)(2)
            << ',' << sys.variable_name(v) << ',' << node->dof_number(sys.number(),v,0) << '\n';
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
  auto jacobian = sys.matrix->zero_clone();
  _problem_ptr->computeResidual(*sys.solution,*residual,0);
  residual->print_matlab(stem + "_residual.m");
  _problem_ptr->computeJacobian(*sys.solution,*jacobian,0);
  jacobian->print_matlab(stem + "_jacobian.m");
  // Re-evaluate the original accepted solution last, with original old histories.
  _problem_ptr->computeResidual(*sys.solution,*residual,0);
  writeHistory("after");
}
