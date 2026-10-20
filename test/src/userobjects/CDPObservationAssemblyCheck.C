#include "CDPObservationAssemblyCheck.h"
#include "CDPObservationProblem.h"
#include "NonlinearSystem.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
registerMooseObject("DamSafetyTestApp", CDPObservationAssemblyCheck);
InputParameters CDPObservationAssemblyCheck::validParams()
{
  auto p=GeneralUserObject::validParams();
  p.set<ExecFlagEnum>("execute_on")=EXEC_INITIAL;
  return p;
}
void CDPObservationAssemblyCheck::execute()
{
  if(_t!=0) mooseError("Development assembly check requires time zero");
  auto * problem=dynamic_cast<CDPObservationProblem *>(&_fe_problem);
  if(!problem) mooseError("Development assembly check requires observation problem");
  auto & nl=_fe_problem.getNonlinearSystem(0);
  auto & system=dynamic_cast<libMesh::NonlinearImplicitSystem &>(nl.system());
  auto matrix=system.get_system_matrix().zero_clone();
  problem->attemptEvent("development_jacobian_before");
  _fe_problem.computeJacobian(*nl.currentSolution(),*matrix,0);
  problem->attemptEvent("development_jacobian_after");
}
