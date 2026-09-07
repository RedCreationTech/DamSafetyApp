#include "CDPTrialProblem.h"
#include "CDPAssemblyProbe.h"
#include "NonlinearSystem.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/numeric_vector.h"
#include <fstream>
#include <iomanip>
registerMooseObject("DamSafetyApp", CDPTrialProblem);
InputParameters CDPTrialProblem::validParams()
{
  auto p = FEProblem::validParams();
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
  if (_observed_snes != snes)
  {
    if (SNESMonitorSet(snes, monitor, this, nullptr)) mooseError("Cannot attach passive SNES monitor");
    _observed_snes = snes;
  }
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
  if (!p.capturing()) return 0;
  // SNES monitor reports the accepted iterate; residual calls alone do not.
  Vec x;
  PetscCall(SNESGetSolution(snes, &x));
  const auto stem=p._prefix+"_accepted_eval"+std::to_string(p._evaluation)+"_it"+std::to_string(iteration);
  PetscViewer viewer;
  PetscCall(PetscViewerASCIIOpen(PetscObjectComm((PetscObject)snes), (stem+"_u.m").c_str(), &viewer));
  PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB));
  PetscCall(VecView(x,viewer));
  PetscCall(PetscViewerDestroy(&viewer));
  if (p.processor_id()==0)
  {
    std::ofstream f(p._prefix+"_accepted.csv",std::ios::app);
    f << std::setprecision(17) << p._evaluation << ',' << p.time() << ',' << iteration << ',' << norm << '\n';
    if (!f) return PETSC_ERR_FILE_WRITE;
  }
  return 0;
}
