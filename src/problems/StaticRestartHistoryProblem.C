#include "StaticRestartHistoryProblem.h"
#include "StaticRestartHistory.h"
#include "AuxiliarySystem.h"
#include "MooseApp.h"
#include "MooseMesh.h"
#include "MooseVariableFieldBase.h"
#include "NonlinearSystemBase.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"
#include <cmath>

registerMooseObject("DamSafetyApp", StaticRestartHistoryProblem);

InputParameters
StaticRestartHistoryProblem::validParams()
{
  auto params = FEProblem::validParams();
  params.addClassDescription("Opt-in static restart kinematic histories; preserves current fields and CDP material storage.");
  params.addParam<bool>("initialize_static_history", false,
                        "After a static checkpoint restart, copy current solution into existing old/older histories once");
  params.addRangeCheckedParam<Real>("static_restart_time", 1.0, "static_restart_time >= 0",
                                   "Exact total time of the accepted static checkpoint");
  return params;
}

StaticRestartHistoryProblem::StaticRestartHistoryProblem(const InputParameters & parameters)
  : FEProblem(parameters),
    _static_history_initialized(declareRestartableData<bool>("static_history_initialized", false))
{
}

void
StaticRestartHistoryProblem::initialSetup()
{
  // The base completes checkpoint loading, ICs, INITIAL objects and deferred restore.
  FEProblem::initialSetup();
  if (!getParam<bool>("initialize_static_history") || _app.isRecovering())
    return;
  if (!_app.isRestarting() || _static_history_initialized ||
      !std::isfinite(time()) || time() != getParam<Real>("static_restart_time"))
    paramError("initialize_static_history", "Requires a not-yet-initialized static checkpoint at static_restart_time; never restart a dynamic state as rest");
  if (numSolverSystems() != 1 || numNonlinearSystems() != 1)
    mooseError("StaticRestartHistoryProblem requires one nonlinear system");

  auto & nl = getNonlinearSystemBase(0);
  auto & aux = getAuxiliarySystem();
  if (!nl.hasSolutionState(1) || !nl.hasSolutionState(2) || !aux.hasSolutionState(1))
    mooseError("Static restart requires the original displacement old/older and auxiliary old histories");

  // These are H39's velocity fields. An IC/configuration error must not be hidden
  // by silently zeroing current fields after the checkpoint has been restored.
  unsigned int invalid_velocity = 0;
  for (const auto * name : {"vel_x", "vel_y"})
  {
    const auto & var = getVariable(0, name);
    if (var.kind() != Moose::VAR_AUXILIARY || &var.sys() != &aux)
      mooseError("Static restart requires auxiliary velocity field ", name);
    for (const auto * node : mesh().getMesh().local_node_ptr_range())
      if (node->n_dofs(aux.system().number(), var.number()))
      {
        const auto dof = node->dof_number(aux.system().number(), var.number(), 0);
        if (dof >= aux.solution().first_local_index() && dof < aux.solution().last_local_index() &&
            aux.solution()(dof) != 0.0)
          invalid_velocity = 1;
      }
  }
  comm().max(invalid_velocity);
  if (invalid_velocity)
    mooseError("Static restart requires zero initial velocity on every MPI rank");

  // Do not call advanceState/copyOldSolutions: they also shift or preserve a ramp.
  StaticRestartHistory::initialize(nl.solution(), nl.solutionOld(), &nl.solutionOlder());
  nl.solutionOld().close();
  nl.solutionOlder().close();
  auto * aux_older = aux.hasSolutionState(2) ? &aux.solutionOlder() : nullptr;
  StaticRestartHistory::initialize(aux.solution(), aux.solutionOld(), aux_older);
  aux.solutionOld().close();
  if (aux_older)
    aux_older->close();
  _static_history_initialized = true;
}
