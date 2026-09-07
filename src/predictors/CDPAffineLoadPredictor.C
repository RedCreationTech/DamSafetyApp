#include "CDPAffineLoadPredictor.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"
#include "MooseMesh.h"
#include "Function.h"
#include "libmesh/system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include <fstream>
#include <iomanip>

registerMooseObject("DamSafetyApp", CDPAffineLoadPredictor);
InputParameters CDPAffineLoadPredictor::validParams()
{
  auto p = Predictor::validParams();
  p.set<Real>("scale") = 1.0;
  p.addParam<bool>("enabled", false, "Explicit opt-in; disabled leaves the initial guess unchanged");
  p.addRequiredParam<std::string>("variable", "First-order nodal Z displacement variable");
  p.addRequiredParam<FunctionName>("top_displacement", "Same time-only displacement function as the top BC");
  p.addRequiredParam<Real>("z_bottom", "Fixed bottom Z coordinate");
  p.addRequiredParam<Real>("z_top", "Prescribed top Z coordinate");
  p.addClassDescription("Add an affine load increment to the previous accepted displacement; no symmetry projection");
  return p;
}
CDPAffineLoadPredictor::CDPAffineLoadPredictor(const InputParameters & p)
  : Predictor(p), _enabled(getParam<bool>("enabled")),
    _variable(getParam<std::string>("variable")),
    _top_displacement(_fe_problem.getFunction(getParam<FunctionName>("top_displacement"))),
    _z_bottom(getParam<Real>("z_bottom")), _z_top(getParam<Real>("z_top"))
{
  if (_z_top <= _z_bottom) paramError("z_top", "Must exceed z_bottom");
  if (_scale != 1.0) paramError("scale", "This boundary-consistent predictor requires scale=1");
}
bool CDPAffineLoadPredictor::shouldApply()
{
  return _enabled && Predictor::shouldApply();
}
void CDPAffineLoadPredictor::apply(NumericVector<Number> & sln)
{
  auto & sys = _nl.system();
  const auto var = sys.variable_number(_variable);
  const auto type = sys.variable_type(var);
  if (type.family != libMesh::LAGRANGE || type.order != libMesh::FIRST)
    mooseError("CDPAffineLoadPredictor requires a FIRST LAGRANGE nodal displacement");
  const Point top(0, 0, _z_top);
  const Real increment = _top_displacement.value(_fe_problem.time(), top) -
                         _top_displacement.value(_fe_problem.timeOld(), top);
  // Always restart from the accepted solution, including retries after cutback.
  // X/Y and every pre-existing asymmetric displacement component are retained.
  sln = _solution_old;
  for (const auto * node : _nl.mesh().getMesh().local_node_ptr_range())
  {
    if (node->n_dofs(_nl.number(), var) != 1)
      mooseError("CDPAffineLoadPredictor requires one selected DOF at every node");
    const Real fraction = ((*node)(2) - _z_bottom) / (_z_top - _z_bottom);
    if (fraction < -1e-10 || fraction > 1.0 + 1e-10)
      mooseError("Node outside the declared predictor height");
    const auto dof = node->dof_number(_nl.number(), var, 0);
    sln.set(dof, _solution_old(dof) + fraction * increment);
  }
  sln.close();
  if (processor_id() == 0)
  {
    std::ofstream f("affine_predictor.csv", std::ios::app);
    f << std::setprecision(17) << _t_step << ',' << _fe_problem.timeOld() << ','
      << _fe_problem.time() << ',' << increment << '\n';
    if (!f) mooseError("Cannot write affine predictor audit");
  }
}
