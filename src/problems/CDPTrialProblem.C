#include "CDPTrialProblem.h"
#include "CDPAssemblyProbe.h"
#include "NonlinearSystem.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "MaterialPropertyStorage.h"
#include "MooseMesh.h"
#include "RankTwoTensor.h"
#include <fstream>
#include <iomanip>
registerMooseObject("DamSafetyApp", CDPTrialProblem);
InputParameters CDPTrialProblem::validParams()
{
  auto p = FEProblem::validParams();
  p.addParam<bool>("audit_directions", false, "One-shot intrusive frozen-state derivative audit; never commits history");
  p.addParam<Real>("audit_time", 0.57, "Exact attempted time for derivative audit");
  p.addParam<unsigned int>("audit_iteration", 1, "SNES accepted iterate for derivative audit");
  p.addParam<bool>("capture_actual_jacobian", false, "Passively write actual assembled solver matrices in the trial window");
  p.addParam<Real>("trial_start", 0.425, "First observed actual residual time");
  p.addParam<Real>("trial_end", 0.430, "Last observed actual residual time");
  p.addParam<std::string>("trial_prefix", "trial", "Output prefix");
  return p;
}
CDPTrialProblem::CDPTrialProblem(const InputParameters & p)
  : FEProblem(p), _start(getParam<Real>("trial_start")), _end(getParam<Real>("trial_end")),
    _prefix(getParam<std::string>("trial_prefix")) {}
bool CDPTrialProblem::capturing() const { return time() >= _start-1e-9 && time() <= _end+1e-9; }
void CDPTrialProblem::computeResidual(const NumericVector<Number> & x,
                                      NumericVector<Number> & r, unsigned int n)
{
  auto * solver = dynamic_cast<libMesh::PetscNonlinearSolver<Number> *>(getNonlinearSystem(n).nonlinearSolver());
  if (!solver) mooseError("CDPTrialProblem requires PETSc");
  SNES snes = solver->snes();
  // SNES objects may be recreated at the same address between time steps.
  // PETSc deduplicates equal monitor/context/destroy triples. Re-registering
  // avoids a rank-dependent stale pointer cache and collective viewer hangs.
  if (SNESMonitorSet(snes, monitor, this, nullptr))
    mooseError("Cannot attach passive SNES monitor");
  if (!capturing()) { FEProblem::computeResidual(x,r,n); return; }
  const auto id = ++_evaluation;
  const auto stem = _prefix + "_eval" + std::to_string(id);
  PetscInt iteration = -1;
  if (SNESGetIterationNumber(snes, &iteration)) mooseError("Cannot query SNES iteration");
  x.print_matlab(stem + "_u.m");
  CDPAssemblyProbe::beginCapture(stem + "_ip_rank" + std::to_string(processor_id()) + ".csv", false);
  try { FEProblem::computeResidual(x,r,n); }
  catch (...) { CDPAssemblyProbe::endCapture(); throw; }
  CDPAssemblyProbe::endCapture();
  r.print_matlab(stem + "_r.m");
  if (processor_id()==0)
  {
    std::ofstream f(_prefix + "_evaluations.csv", std::ios::app);
    f << std::setprecision(17) << id << ',' << time() << ',' << dt() << ',' << iteration << '\n';
    if (!f) mooseError("Cannot write trial index");
  }
}
PetscErrorCode CDPTrialProblem::monitor(SNES snes, PetscInt iteration, PetscReal norm, void * context)
{
  auto & p = *static_cast<CDPTrialProblem *>(context);
  if (!p.capturing()) return PETSC_SUCCESS;
  // SNES monitor reports the accepted iterate; residual calls alone do not.
  Vec x;
  PetscCall(SNESGetSolution(snes, &x));
  const auto stem=p._prefix+"_accepted_eval"+std::to_string(p._evaluation)+"_it"+std::to_string(iteration);
  PetscViewer viewer;
  PetscCall(PetscViewerASCIIOpen(PetscObjectComm((PetscObject)snes), (stem+"_u.m").c_str(), &viewer));
  PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB));
  PetscCall(VecView(x,viewer));
  PetscCall(PetscViewerDestroy(&viewer));
  if (iteration > 0 && p.getParam<bool>("capture_actual_jacobian"))
  {
    // Observe the operator actually used for the completed linear solve.
    // Never print/close the live matrix inside FEProblem assembly.
    KSP ksp;
    Mat op, pc, copy;
    PetscCall(SNESGetKSP(snes, &ksp));
    PetscCall(KSPGetOperators(ksp, &op, &pc));
    PetscCall(MatDuplicate(pc, MAT_COPY_VALUES, &copy));
    const auto matrix_stem = p._prefix + "_actual_jac_eval" + std::to_string(p._evaluation);
    PetscCall(PetscViewerASCIIOpen(PetscObjectComm((PetscObject)snes), (matrix_stem+"_J.m").c_str(), &viewer));
    PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(MatView(copy, viewer));
    PetscCall(PetscViewerDestroy(&viewer));
    PetscCall(MatDestroy(&copy));
    if (p.processor_id()==0)
    {
      std::ofstream f(p._prefix+"_actual_jacobians.csv",std::ios::app);
      f << std::setprecision(17) << p._evaluation << ',' << p.time() << ',' << p.dt() << ',' << iteration << '\n';
      if (!f) return PETSC_ERR_FILE_WRITE;
    }
  }
  if (p.processor_id()==0)
  {
    std::ofstream f(p._prefix+"_accepted.csv",std::ios::app);
    f << std::setprecision(17) << p._evaluation << ',' << p.time() << ',' << iteration << ',' << norm << '\n';
    if (!f) return PETSC_ERR_FILE_WRITE;
  }
  if (p.getParam<bool>("audit_directions") && !p._directions_done &&
      std::abs(p.time()-p.getParam<Real>("audit_time")) < 1e-9 &&
      iteration == static_cast<PetscInt>(p.getParam<unsigned int>("audit_iteration")))
  {
    p._directions_done = true;
    p.auditDirections(snes);
  }
  return PETSC_SUCCESS;
}

void CDPTrialProblem::auditDirections(SNES snes)
{
  auto & nl = getNonlinearSystem(0);
  auto & sys = dynamic_cast<libMesh::NonlinearImplicitSystem &>(nl.system());
  const auto * original = nl.currentSolution();
  auto x = original->clone();
  auto residual = x->zero_clone();
  auto jacobian = sys.get_system_matrix().zero_clone();
  const std::string stem = _prefix + "_direction";
  const std::string rank = "_rank" + std::to_string(processor_id()) + ".csv";
  auto history = [&](const std::string & label) {
    const auto & storage = getMaterialPropertyStorage();
    std::ofstream f(stem+"_history_"+label+rank);
    f << std::setprecision(17) << "element,side,state,property,qp,component,value\n";
    for (auto state : storage.stateIndexRange())
      for (const auto & entry : storage.props(state))
        for (const auto & side : entry.second)
          for (unsigned int id=0; id<side.second.size(); ++id)
          {
            const auto * value=side.second.queryValue(id);
            if (!value) continue;
            const auto name=storage.queryStatefulPropName(storage.statefulProps().at(id));
            const auto * scalar=dynamic_cast<const MaterialProperty<Real> *>(value);
            const auto * tensor=dynamic_cast<const MaterialProperty<RankTwoTensor> *>(value);
            if (!name || (!scalar && !tensor)) mooseError("Unsupported audit history");
            for (unsigned int qp=0;qp<value->size();++qp)
              for (unsigned int c=0;c<(scalar?1u:9u);++c)
                f << entry.first->id() << ',' << side.first << ',' << state << ',' << *name << ',' << qp << ',' << c << ',' << (scalar?(*scalar)[qp]:(*tensor)[qp](c/3,c%3)) << '\n';
          }
    if (!f) mooseError("Cannot write direction history");
  };
  auto evaluate = [&](const NumericVector<Number> & u, const std::string & label) {
    CDPAssemblyProbe::beginCapture(stem+"_"+label+rank,false);
    try { FEProblem::computeResidual(u,*residual,0); }
    catch (...) { CDPAssemblyProbe::endCapture(); throw; }
    CDPAssemblyProbe::endCapture();
    residual->print_matlab(stem+"_"+label+"_r.m");
  };
  x->print_matlab(stem+"_u.m");
  nl.solutionOld().print_matlab(stem+"_old_u.m");
  history("before");
  evaluate(*x,"base0");
  evaluate(*x,"base1");
  CDPAssemblyProbe::beginCapture(stem+"_jacobian"+rank,true);
  try { FEProblem::computeJacobian(*x,*jacobian,0); }
  catch (...) { CDPAssemblyProbe::endCapture(); throw; }
  CDPAssemblyProbe::endCapture();
  jacobian->print_matlab(stem+"_J.m");
  evaluate(*x,"base2");
  auto direction=x->zero_clone();
  Vec update;
  if (SNESGetSolutionUpdate(snes,&update)) mooseError("Missing Newton update");
  const PetscScalar * a;
  PetscInt lo,hi;
  if (VecGetOwnershipRange(update,&lo,&hi) || VecGetArrayRead(update,&a)) mooseError("Cannot read Newton update");
  for (PetscInt i=lo;i<hi;++i) direction->set(i,a[i-lo]);
  if (VecRestoreArrayRead(update,&a)) mooseError("Cannot restore Newton update view");
  direction->close();
  // Boundary-compatible Newton direction; the prescribed planes are fixed at this time.
  std::vector<Number> all;
  direction->localize(all);
  const auto & mesh = this->mesh().getMesh();
  Real zmin=1e30,zmax=-1e30;
  for (const auto * node : mesh.local_node_ptr_range()) { zmin=std::min(zmin,(*node)(2));zmax=std::max(zmax,(*node)(2)); }
  comm().min(zmin);comm().max(zmax);
  for (const auto * node : mesh.local_node_ptr_range())
    if (std::abs((*node)(2)-zmin)<1e-12 || std::abs((*node)(2)-zmax)<1e-12)
      for (unsigned int v=0;v<sys.n_vars();++v)
        if (node->n_dofs(sys.number(),v)) direction->set(node->dof_number(sys.number(),v,0),0.);
  direction->close();
  for (unsigned int kind=0;kind<2;++kind)
  {
    if (kind==1)
    {
      // Smooth Z-reflection-odd vector: Z component even about midplane,
      // zero on both prescribed planes. It does not alter boundary data.
      direction->zero();
      const auto v=sys.variable_number("disp_z");
      for (const auto * node : mesh.local_node_ptr_range())
      {
        const Real z=((*node)(2)-zmin)/(zmax-zmin);
        direction->set(node->dof_number(sys.number(),v,0),4*z*(1-z));
      }
      direction->close();
    }
    const auto scale=direction->linfty_norm();
    if (scale==0.) mooseError("Zero audit direction");
    direction->scale(1./scale);
    direction->close();
    const auto label=std::string(kind==0?"newton":"odd");
    direction->print_matlab(stem+"_"+label+"_v.m");
    auto jv=x->zero_clone();jacobian->vector_mult(*jv,*direction);
    jv->print_matlab(stem+"_"+label+"_Jv.m");
    unsigned int k=0;
    for (const Real h : {1e-8,1e-9,1e-10,1e-11})
    {
      for (const Real sign : {1.,-1.})
      {
        auto trial=x->clone();trial->add(sign*h,*direction);trial->close();
        evaluate(*trial,label+"_h"+std::to_string(k)+(sign>0?"_plus":"_minus"));
        // Restore the current solution pointer before each scratch vector dies.
        evaluate(*x,label+"_h"+std::to_string(k)+(sign>0?"_resetp":"_resetm"));
      }
      ++k;
    }
  }
  evaluate(*original,"restored");
  history("after");
  if (processor_id()==0)
  {
    std::ofstream f(stem+"_identity.csv");
    f << "time,dt,evaluation,h0,h1,h2,h3\n" << std::setprecision(17) << time() << ',' << dt() << ',' << _evaluation << ",1e-8,1e-9,1e-10,1e-11\n";
  }
}
