#include "CDPObservationProblem.h"

#include "CDPObservationContextGuard.h"
#include "CDPObservationFingerprint.h"
#include "Assembly.h"
#include "KernelBase.h"
#include "AuxiliarySystem.h"
#include "MaterialData.h"
#include "MaterialPropertyStorage.h"
#include "MooseMesh.h"
#include "MooseMeshUtils.h"
#include "NonlinearSystem.h"
#include "RankTwoTensor.h"
#include "libmesh/elem.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_vector.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>


namespace
{
std::map<std::string, std::vector<std::string>>
parseGroupTerms(const std::vector<std::string> & spec)
{
  std::map<std::string, std::vector<std::string>> terms;
  for (const auto & entry : spec)
  {
    const auto colon = entry.find(':');
    if (colon == std::string::npos || colon == 0)
      mooseError("observation_residual_group_terms entries must look like 'term:name1,name2': ",
                 entry);
    std::vector<std::string> names;
    std::stringstream rest(entry.substr(colon + 1));
    std::string name;
    while (std::getline(rest, name, ','))
      if (!name.empty())
        names.push_back(name);
    if (names.empty())
      mooseError("observation_residual_group_terms term has no kernels: ", entry);
    terms[entry.substr(0, colon)] = names;
  }
  return terms;
}
} // namespace

registerMooseObject("DamSafetyApp", CDPObservationProblem);

namespace
{
std::string
quoted(const std::string & value)
{
  std::string result = "\"";
  for (const char c : value)
  {
    if (c == '"')
      result += '"';
    result += c;
  }
  return result + '"';
}
}

InputParameters
CDPObservationProblem::validParams()
{
  auto parameters = FEProblem::validParams();
  parameters.addParam<bool>(
      "observation_enabled", false, "Enable passive recording without additional evaluations");
  parameters.addParam<bool>(
      "trace_state_swap", false, "Record state and completed candidates around existing swap-back");
  parameters.addParam<Real>("observation_start", 3.018, "First observed cumulative time");
  parameters.addParam<Real>("observation_end", 3.202, "Last observed cumulative time");
  parameters.addParam<Real>("earthquake_offset", 1.0, "Cumulative time origin of EQ");
  parameters.addParam<std::string>(
      "observation_prefix", "ds04g-observation", "Per-rank output file prefix");
  parameters.addParam<unsigned long long>("observation_byte_limit",
                                          65536000,
                                          "Per-rank detailed-output byte cap");
  parameters.addParam<std::vector<unsigned int>>(
      "observation_labels", {}, "Source element labels");
  parameters.addParam<std::vector<Real>>(
      "observation_vertices", {}, "Four x/y vertices per label in SI");
  parameters.addParam<bool>(
      "observation_residual_groups", false, "Record grouped nodal residuals at accepted steps in window");
  parameters.addParam<std::vector<std::string>>(
      "observation_residual_group_terms",
      std::vector<std::string>{"entity:TM_all0,TM_all1",
                               "mass_base:TM_all_inertia_x,TM_all_inertia_y",
                               "mass_correction:ds04g_sp1_mass_0,ds04g_sp1_mass_1",
                               "gravity_base:gravity_y",
                               "gravity_correction:ds04g_sp1_gravity_1",
                               "hg:ds04g_sp1_hourglass_0,ds04g_sp1_hourglass_1",
                               "entity_rate:ds04g_all_velocity_damping_x,ds04g_all_velocity_damping_y"},
      "Term to kernel names as 'term:name1,name2'; derived union terms (mass/gravity/all) are built automatically");
  parameters.addParam<std::vector<unsigned int>>(
      "observation_residual_nodes", std::vector<unsigned int>{35, 167}, "Extra node ids for grouped residuals");
  parameters.addParam<BoundaryName>(
      "observation_residual_boundary", "DAM_BASE__Dam_1", "Boundary whose nodes give the net-force sums");
  return parameters;
}

CDPObservationProblem::CDPObservationProblem(const InputParameters & parameters)
  : FEProblem(parameters),
    _enabled(getParam<bool>("observation_enabled")),
    _trace_swap(getParam<bool>("trace_state_swap")),
    _start(getParam<Real>("observation_start")),
    _end(getParam<Real>("observation_end")),
    _prefix(getParam<std::string>("observation_prefix")),
    _limit(getParam<unsigned long long>("observation_byte_limit")),
    _res_groups(getParam<bool>("observation_residual_groups")),
    _res_group_terms(parseGroupTerms(getParam<std::vector<std::string>>(
        "observation_residual_group_terms"))),
    _res_extra_nodes(getParam<std::vector<unsigned int>>("observation_residual_nodes")),
    _res_boundary(getParam<BoundaryName>("observation_residual_boundary"))
{
  if (_enabled && libMesh::n_threads() != 1)
    paramError("observation_enabled", "Requires one thread per MPI rank");
  if (_start > _end)
    paramError("observation_start", "Must not exceed observation_end");
}

CDPObservationProblem::~CDPObservationProblem()
{
  if (_enabled)
    status();
  if (CDPObservation::configured == this)
    CDPObservation::configured = nullptr;
}

void
CDPObservationProblem::initialSetup()
{
  if (!_enabled)
  {
    FEProblem::initialSetup();
    return;
  }

  const auto & labels = getParam<std::vector<unsigned int>>("observation_labels");
  const auto & vertices = getParam<std::vector<Real>>("observation_vertices");
  if (labels.empty() || vertices.size() != labels.size() * 8)
    paramError("observation_vertices", "Expected eight coordinates per label");
  if (std::set<unsigned int>(labels.begin(), labels.end()).size() != labels.size())
    paramError("observation_labels", "Labels must be unique");
  for (const auto value : vertices)
    if (!std::isfinite(value))
      paramError("observation_vertices", "Coordinates must be finite");

  for (unsigned int k = 0; k < labels.size(); ++k)
  {
    unsigned int matches = 0;
    for (const auto * element : mesh().getMesh().active_local_element_ptr_range())
    {
      if (element->dim() != 2 || element->n_nodes() != 4)
        continue;
      unsigned int hits = 0;
      bool used[4] = {false, false, false, false};
      bool ambiguous = false;
      for (unsigned int i = 0; i < 4; ++i)
      {
        unsigned int count = 0;
        unsigned int match = 0;
        for (unsigned int j = 0; j < 4; ++j)
          if (std::abs(element->point(i)(0) - vertices[k * 8 + j * 2]) < 1e-7 &&
              std::abs(element->point(i)(1) - vertices[k * 8 + j * 2 + 1]) < 1e-7)
          {
            ++count;
            match = j;
          }
        if (count != 1 || used[match])
        {
          ambiguous = true;
          break;
        }
        used[match] = true;
        ++hits;
      }
      if (!ambiguous && hits == 4)
      {
        ++matches;
        _targets.emplace(element->id(), labels[k]);
      }
    }
    comm().sum(matches);
    if (matches != 1)
      paramError("observation_vertices", "Source element mapping must be unique: ", labels[k]);
  }

  if (CDPObservation::configured)
    mooseError("Only one observation problem per process is supported");
  CDPObservation::configured = this;
  FEProblem::initialSetup();
  if (_res_groups)
  {
    std::set<dof_id_type> nodes;
    for (const auto & target : _targets)
    {
      const auto * element = mesh().getMesh().elem_ptr(target.first);
      for (unsigned int i = 0; i < element->n_nodes(); ++i)
        nodes.insert(element->node_id(i));
    }
    for (const auto id : _res_extra_nodes)
    {
      if (!mesh().queryNodePtr(id))
        paramError("observation_residual_nodes", "Node not found in mesh: ", id);
      nodes.insert(id);
    }
    const auto bid = MooseMeshUtils::getBoundaryID(_res_boundary, mesh().getMesh());
    if (bid == Moose::INVALID_BOUNDARY_ID)
      paramError("observation_residual_boundary", "Boundary not found: ", _res_boundary);
    const auto & binfo = mesh().getMesh().get_boundary_info();
    for (const auto * node : mesh().getMesh().node_ptr_range())
    {
      std::vector<boundary_id_type> ids;
      binfo.boundary_ids(node, ids);
      if (std::find(ids.begin(), ids.end(), bid) != ids.end())
      {
        _res_base_nodes.insert(node->id());
        nodes.insert(node->id());
      }
    }
    _res_target_nodes.assign(nodes.begin(), nodes.end());
    auto & nl = getNonlinearSystem(0);
    std::set<std::string> existing;
    for (const auto & k : nl.getKernelWarehouse().getObjects())
      existing.insert(k->name());
    for (const auto & term : _res_group_terms)
      for (const auto & name : term.second)
        if (!existing.count(name))
          paramError("observation_residual_group_terms",
                     "Kernel not found for term ", term.first, ": ", name);
  }
  for (const auto & target : _targets)
  {
    _element = target.first;
    const auto * element = mesh().getMesh().elem_ptr(target.first);
    for (unsigned int i = 0; i < 4; ++i)
      write("probe-map",
            "vertex",
            "node_id_x_y_z",
            {static_cast<double>(element->node_id(i)),
             element->point(i)(0),
             element->point(i)(1),
             element->point(i)(2)},
            {});
  }
  _element = 0;
  attemptEvent("initial_state");
  status();
}

bool
CDPObservationProblem::window() const
{
  return _enabled && time() >= _start - 1e-9 && time() <= _end + 1e-9;
}

void
CDPObservationProblem::setExecFlag(const ExecFlagType & type) noexcept
{
  try
  {
    _exec_flag = type.name();
  }
  catch (...)
  {
    _io_failed = true;
  }
}

void
CDPObservationProblem::execute(const ExecFlagType & exec_type)
{
  if (!_enabled)
  {
    FEProblem::execute(exec_type);
    return;
  }
  CDPObservationFEContextGuard restore(_evaluation, _context, _global_iteration);
  CDPObservationContextGuard restore_point(_element, _qp, _call);
  ++_execution;
  _evaluation = 0;
  _context = "execute";
  _global_iteration = -1;
  _element = 0;
  _qp = 0;
  _call = 0;
  setExecFlag(exec_type);
  write("attempt-events", "execute_begin", "exec_flag", {}, _exec_flag);
  try
  {
    FEProblem::execute(exec_type);
  }
  catch (...)
  {
    write("attempt-events", "execute_exception", "exec_flag", {}, _exec_flag);
    if (window())
      status();
    throw;
  }
  write("attempt-events", "execute_complete", "exec_flag", {}, _exec_flag);
  if (window())
    status();
}

void
CDPObservationProblem::outputStep(ExecFlagType type)
{
  if (!_enabled)
  {
    FEProblem::outputStep(type);
    return;
  }
  CDPObservationFEContextGuard restore(_evaluation, _context, _global_iteration);
  CDPObservationContextGuard restore_point(_element, _qp, _call);
  ++_output;
  _evaluation = 0;
  _context = "output";
  _global_iteration = -1;
  _element = 0;
  _qp = 0;
  _call = 0;
  setExecFlag(type);
  write("attempt-events", "output_begin", "output_id", {static_cast<double>(_output)}, _exec_flag);
  try
  {
    FEProblem::outputStep(type);
  }
  catch (...)
  {
    write("attempt-events", "output_exception", "output_id", {static_cast<double>(_output)}, _exec_flag);
    if (window())
      status();
    throw;
  }
  write("attempt-events", "output_complete", "output_id", {static_cast<double>(_output)}, _exec_flag);
  if (window())
    status();
}

void
CDPObservationProblem::advanceState()
{
  if (!_enabled)
  {
    FEProblem::advanceState();
    return;
  }
  if (window())
    attemptEvent("state_advance_before");
  FEProblem::advanceState();
  if (window())
    attemptEvent("state_advance_after");
}

void
CDPObservationProblem::materialCompleted() noexcept
{
  if (!_enabled || !_trace_swap || !window())
    return;
  try
  {
    if (!_targets.count(_element) || _qp != 0 || !_call)
    {
      _swap_incomplete = true;
      return;
    }
    _swap_tracker.completed(_element, _attempt, _execution, _evaluation, _call);
  }
  catch (...)
  {
    _swap_incomplete = true;
  }
}

void
CDPObservationProblem::swapBackMaterials(const THREAD_ID tid)
{
  if (!_enabled || !_trace_swap || !window())
  {
    FEProblem::swapBackMaterials(tid);
    return;
  }

  const Elem * element = nullptr;
  bool was_swapped = false;
  try
  {
    if (tid != 0)
      _swap_incomplete = true;
    else
    {
      element = assembly(tid, 0).elem();
      was_swapped = getMaterialPropertyStorage().getMaterialData(tid).isSwapped();
    }
  }
  catch (...)
  {
    _swap_incomplete = true;
  }

  FEProblem::swapBackMaterials(tid);
  if (!element || !_targets.count(element->id()))
    return;

  try
  {
    const auto batch = _swap_tracker.take(element->id());
    ++_swap_events;
    if (!was_swapped || batch.overflow || batch.calls.size() != 1)
      _swap_incomplete = true;
    if (batch.calls.empty())
      ++_swap_empty;
    CDPObservationContextGuard restore_point(_element, _qp, _call);
    CDPObservationFEContextGuard restore(_evaluation, _context, _global_iteration);
    _element = element->id();
    _qp = 0;
    _call = 0;
    _evaluation = 0;
    _context = "swapback";
    _global_iteration = -1;
    write("swapback-events",
          "swapback_after",
          "was_swapped_count_overflow",
          {was_swapped ? 1.0 : 0.0,
           static_cast<double>(batch.calls.size()),
           batch.overflow ? 1.0 : 0.0},
          "completed candidates since previous swapback; not an accepted-call verdict");
    for (const auto & candidate : batch.calls)
      write("swapback-events",
            "swapback_candidate",
            "attempt_execution_evaluation_call",
            {static_cast<double>(candidate.attempt),
             static_cast<double>(candidate.execution),
             static_cast<double>(candidate.evaluation),
             static_cast<double>(candidate.call)},
            "candidate must be checked against stored state and subsequent advance");
    history("material_swapback_after");
  }
  catch (...)
  {
    _swap_incomplete = true;
    _history_incomplete = true;
  }
}

bool
CDPObservationProblem::enterPoint(std::uint64_t element, unsigned int qp, bool new_call) noexcept
{
  if (!_enabled)
    return false;
  try
  {
    if (new_call)
      _call = ++_next_call;
    else
      _call = 0;
    if (!window() || !_targets.count(element) || qp != 0)
      return false;
    _element = element;
    _qp = qp;
    if (new_call)
      _last_call[{element, qp}] = {_attempt, _execution, _evaluation, _call};
    else
    {
      const auto it = _last_call.find({element, qp});
      if (it != _last_call.end() && it->second.evaluation == _evaluation &&
          it->second.execution == _execution)
        _call = it->second.call;
    }
    return true;
  }
  catch (...)
  {
    _io_failed = true;
    return false;
  }
}

void
CDPObservationProblem::write(const char * stream,
                             const char * event,
                             const char * name,
                             const std::vector<double> & values,
                             const std::string & note) noexcept
{
  try
  {
    if (!_enabled || _io_failed)
      return;
    const auto write_start = std::chrono::steady_clock::now();
    if (_truncated)
      return;
    auto & file = _streams[stream];
    const std::string header =
        file.is_open()
            ? ""
            : "rank,record,time,dt,dt_old,step,attempt,evaluation,execution,output,call,"
              "global_iteration,context,exec_flag,element,source_element,qp,partition,N,substep,"
              "local_solve,local_iteration,event,field,values,note\n";
    std::ostringstream encoded_values;
    encoded_values << std::setprecision(17);
    for (unsigned int i = 0; i < values.size(); ++i)
    {
      if (i)
        encoded_values << ' ';
      encoded_values << values[i];
    }
    const auto & context = CDPObservation::context;
    std::ostringstream row;
    row << std::setprecision(17) << processor_id() << ',' << _record + 1 << ',' << time() << ','
        << dt() << ',' << dtOld() << ',' << timeStep() << ',' << _attempt << ',' << _evaluation
        << ',' << _execution << ',' << _output << ',' << _call << ',' << _global_iteration << ','
        << _context << ',' << _exec_flag << ',' << _element << ','
        << (_targets.count(_element) ? _targets.at(_element) : 0) << ',' << _qp << ','
        << context.partition << ',' << context.N << ',' << context.substep << ','
        << context.local_solve << ',' << context.iteration << ',' << event << ',' << quoted(name)
        << ',' << quoted(encoded_values.str()) << ',' << quoted(note) << '\n';
    const auto text = row.str();
    if (_bytes > _limit || header.size() > _limit - _bytes ||
        text.size() > _limit - _bytes - header.size())
    {
      _truncated = true;
      return;
    }
    if (!file.is_open())
      file.open(_prefix + "-" + stream + "-rank" + std::to_string(processor_id()) + ".csv",
                std::ios::out);
    file << header << text;
    _bytes += header.size() + text.size();
    ++_record;
    ++_counts[stream];
    if (!file)
      _io_failed = true;
    if (std::string(event).find("exception") != std::string::npos)
      file.flush();
    _write_microseconds +=
        std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - write_start)
            .count();
  }
  catch (...)
  {
    _io_failed = true;
  }
}

void
CDPObservationProblem::status() noexcept
{
  try
  {
    for (auto & stream : _streams)
    {
      stream.second.flush();
      if (!stream.second)
        _io_failed = true;
    }
    std::ofstream file(_prefix + "-status-rank" + std::to_string(processor_id()) + ".json");
    file << "{\"schema\":\"ds04g-observation-v2\",\"rank\":" << processor_id()
         << ",\"records\":" << _record << ",\"bytes\":" << _bytes
         << ",\"capture_truncated\":" << (_truncated ? "true" : "false")
         << ",\"io_failed\":" << (_io_failed ? "true" : "false")
         << ",\"history_incomplete\":" << (_history_incomplete ? "true" : "false")
         << ",\"target_elements_on_rank\":" << _targets.size()
         << ",\"write_microseconds\":" << _write_microseconds
         << ",\"swap_trace_enabled\":" << (_trace_swap ? "true" : "false")
         << ",\"swap_trace_events\":" << _swap_events << ",\"swap_trace_empty\":" << _swap_empty
         << ",\"swap_trace_incomplete\":" << (_swap_incomplete ? "true" : "false")
         << ",\"residual_groups_enabled\":" << (_res_groups ? "true" : "false")
         << ",\"residual_groups_rows\":" << _res_rows
         << ",\"residual_groups_contaminated\":" << (_res_contaminated ? "true" : "false")
         << ",\"executions\":" << _execution << ",\"outputs\":" << _output
         << ",\"material_calls\":" << _next_call << ",\"streams\":{";
    bool comma = false;
    for (const auto & item : _counts)
    {
      if (comma)
        file << ',';
      comma = true;
      file << quoted(item.first) << ':' << item.second;
    }
    file << "}}\n";
  }
  catch (...)
  {
    _io_failed = true;
  }
}

void
CDPObservationProblem::history(const char * event)
{
  const auto & storage = getMaterialPropertyStorage();
  for (const auto state : storage.stateIndexRange())
    for (const auto & entry : storage.props(state))
    {
      if (!_targets.count(entry.first->id()))
        continue;
      _element = entry.first->id();
      for (const auto & side : entry.second)
        for (unsigned int id = 0; id < side.second.size(); ++id)
        {
          const auto * value = side.second.queryValue(id);
          if (!value)
            continue;
          const auto name = storage.queryStatefulPropName(storage.statefulProps().at(id));
          const auto * scalar = dynamic_cast<const MaterialProperty<Real> *>(value);
          const auto * tensor = dynamic_cast<const MaterialProperty<RankTwoTensor> *>(value);
          if (!name || (!scalar && !tensor))
          {
            _history_incomplete = true;
            continue;
          }
          if (value->size() == 0)
            continue;
          _qp = 0;
          std::vector<double> values;
          if (scalar)
            values.push_back((*scalar)[0]);
          else
            for (unsigned int component = 0; component < 9; ++component)
              values.push_back((*tensor)[0](component / 3, component % 3));
          const auto field = "state" + std::to_string(state) + "_side" +
                             std::to_string(side.first) + "_" + *name;
          write("history-events",
                event,
                field.c_str(),
                values,
                state == 0 ? "stored_current_unverified" : "stored_accepted_history");
        }
  }
  _element = 0;
  _qp = 0;

  for (auto * system : {static_cast<SystemBase *>(&getNonlinearSystem(0)),
                        static_cast<SystemBase *>(&getAuxiliarySystem())})
  {
    const auto & sys = system->system();
    const auto * current = system->currentSolution();
    const SystemBase & stored = *system;
    const auto * old = stored.hasSolutionState(1) ? &stored.solutionState(1) : nullptr;
    const auto * older = stored.hasSolutionState(2) ? &stored.solutionState(2) : nullptr;
    for (const auto * node : mesh().getMesh().local_node_ptr_range())
    {
      bool selected = false;
      const auto & vertices = getParam<std::vector<Real>>("observation_vertices");
      for (unsigned int j = 0; j < vertices.size(); j += 2)
        if (std::abs((*node)(0) - vertices[j]) < 1e-7 &&
            std::abs((*node)(1) - vertices[j + 1]) < 1e-7)
          selected = true;
      if (!selected)
        continue;
      for (unsigned int variable = 0; variable < sys.n_vars(); ++variable)
      {
        const auto name = sys.variable_name(variable);
        if (name != "disp_x" && name != "disp_y" && name != "vel_x" && name != "vel_y" &&
            name != "accel_x" && name != "accel_y")
          continue;
        if (!node->n_dofs(sys.number(), variable))
          continue;
        const auto dof = node->dof_number(sys.number(), variable, 0);
        if (dof < current->first_local_index() || dof >= current->last_local_index())
          continue;
        write("history-events",
              event,
              name.c_str(),
              {static_cast<double>(node->id()),
               static_cast<double>(dof),
               (*current)(dof),
               old ? (*old)(dof) : std::numeric_limits<double>::quiet_NaN(),
               older ? (*older)(dof) : std::numeric_limits<double>::quiet_NaN()},
              std::string("node_dof_current_old_older; ") +
                  (old ? "old_allocated; " : "old_unallocated; ") +
                  (older ? "older_allocated" : "older_unallocated"));
      }
    }
  }
  _element = 0;
  _qp = 0;
}

std::vector<double>
CDPObservationProblem::acceptedFingerprint()
{
  std::vector<double> values;
  const auto & storage = getMaterialPropertyStorage();
  for (auto state : storage.stateIndexRange())
  {
    if (!state)
      continue; // state 0 is transient and legitimately recomputed by diagnostics
    for (const auto & entry : storage.props(state))
      for (const auto & side : entry.second)
      {
        const auto & props = storage.props(entry.first, side.first, state);
        for (unsigned int id = 0; id < props.size(); ++id)
        {
          const auto * value = props.queryValue(id);
          if (!value)
            continue;
          const auto * scalar = dynamic_cast<const MaterialProperty<Real> *>(value);
          const auto * tensor = dynamic_cast<const MaterialProperty<RankTwoTensor> *>(value);
          if (!scalar && !tensor)
            continue;
          for (unsigned int qp = 0; qp < value->size(); ++qp)
            if (scalar)
              values.push_back((*scalar)[qp]);
            else
              for (unsigned int c = 0; c < 9; ++c)
                values.push_back((*tensor)[qp](c / 3, c % 3));
        }
      }
  }
  for (auto * system : {static_cast<const SystemBase *>(&getNonlinearSystem(0)),
                        static_cast<const SystemBase *>(&getAuxiliarySystem())})
    for (unsigned int state = 1; state < 3; ++state)
    {
      std::vector<Number> data;
      system->solutionState(state).localize(data);
      for (const auto & v : data)
        values.push_back(v);
    }
  values.push_back(dtOld());
  return values;
}

void
CDPObservationProblem::restoreKernelEnable(const std::map<std::string, bool> & enabled) noexcept
{
  try
  {
    for (const auto & k : getNonlinearSystem(0).getKernelWarehouse().getObjects())
      const_cast<bool &>(k->getParam<bool>("enable")) = enabled.at(k->name());
    updateActiveObjects();
  }
  catch (...)
  {
    _io_failed = true;
  }
}

void
CDPObservationProblem::residualGroups(const char * event) noexcept
{
  if (_res_contaminated || _truncated || _io_failed)
    return;
  std::map<std::string, bool> enabled;
  try
  {
    auto & nl = getNonlinearSystem(0);
    const auto & sys = nl.system();
    const auto vdx = sys.variable_number("disp_x");
    const auto vdy = sys.variable_number("disp_y");
    for (const auto & k : nl.getKernelWarehouse().getObjects())
      enabled[k->name()] = k->enabled();
    std::map<std::string, std::vector<std::string>> terms(_res_group_terms.begin(),
                                                          _res_group_terms.end());
    auto unite = [&](const std::string & out, const std::vector<std::string> & parts)
    {
      std::set<std::string> names;
      for (const auto & part : parts)
        if (terms.count(part))
          for (const auto & n : terms.at(part))
            names.insert(n);
      terms[out] = std::vector<std::string>(names.begin(), names.end());
    };
    unite("mass", {"mass_base", "mass_correction"});
    unite("gravity", {"gravity_base", "gravity_correction"});
    unite("all", {"entity", "mass", "gravity", "hg", "entity_rate"});

    const auto before = acceptedFingerprint();
    auto residual = nl.currentSolution()->zero_clone();
    std::vector<dof_id_type> dofs_x, dofs_y;
    std::vector<char> has;
    for (const auto id : _res_target_nodes)
    {
      const auto * node = mesh().queryNodePtr(id);
      const bool hx = node && node->n_dofs(sys.number(), vdx);
      const bool hy = node && node->n_dofs(sys.number(), vdy);
      has.push_back(static_cast<char>((hx ? 1 : 0) | (hy ? 2 : 0)));
      dofs_x.push_back(hx ? node->dof_number(sys.number(), vdx, 0) : DofObject::invalid_id);
      dofs_y.push_back(hy ? node->dof_number(sys.number(), vdy, 0) : DofObject::invalid_id);
    }
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto & term : terms)
    {
      for (const auto & k : nl.getKernelWarehouse().getObjects())
        const_cast<bool &>(k->getParam<bool>("enable")) =
            std::find(term.second.begin(), term.second.end(), k->name()) != term.second.end();
      updateActiveObjects();
      FEProblem::computeResidual(*nl.currentSolution(), *residual, 0);
      std::vector<double> values;
      values.reserve(2 * _res_target_nodes.size() + 2);
      double base_x = 0.0, base_y = 0.0;
      for (std::size_t i = 0; i < _res_target_nodes.size(); ++i)
      {
        const auto id = _res_target_nodes[i];
        const double rx = (has[i] & 1) ? (*residual)(dofs_x[i]) : nan;
        const double ry = (has[i] & 2) ? (*residual)(dofs_y[i]) : nan;
        values.push_back(rx);
        values.push_back(ry);
        if (_res_base_nodes.count(id))
        {
          if (has[i] & 1)
            base_x += rx;
          if (has[i] & 2)
            base_y += ry;
        }
      }
      values.push_back(base_x);
      values.push_back(base_y);
      write("residual-groups",
            event,
            term.first.c_str(),
            values,
            std::string("R_at_accepted_state;nodes=") + std::to_string(_res_target_nodes.size()) +
                ";base_nodes=" + std::to_string(_res_base_nodes.size()) +
                ";value_order=per_target_node_rx_ry_then_net_base_Fx_Fy");
      ++_res_rows;
    }
    restoreKernelEnable(enabled);
    if (acceptedFingerprint() != before)
    {
      _res_contaminated = true;
      write("residual-groups", event, "freeze_check", {1.0},
            "accepted fingerprint changed by grouped residual capture; capture stopped");
    }
  }
  catch (...)
  {
    _res_contaminated = true;
    restoreKernelEnable(enabled);
    write("residual-groups", event, "capture_exception", {1.0},
          "grouped residual capture aborted by exception; capture stopped; solve state untouched");
  }
}

void
CDPObservationProblem::attemptEvent(const char * event) noexcept
{
  if (!_enabled)
    return;
  try
  {
    const std::string event_name(event);
    if (event_name == "attempt_enter")
    {
      ++_attempt;
      _global_iteration = -1;
    }
    if (event_name == "accepted_end" || event_name == "stepper_accept_after")
      _context = "accept";
    else if (event_name.find("state_advance_") == 0)
      _context = "advance";
    else if (event_name.find("reject_") == 0)
      _context = "reject";
    else
      _context = "event";
    _evaluation = 0;
    _global_iteration = -1;
    _element = 0;
    _qp = 0;
    _call = 0;
    CDPObservation::context = {};
    write("attempt-events", event, "time", {time()}, {});
    if (window())
      history(event);
    if (_res_groups && window() && event_name == "state_advance_before")
      residualGroups(event);
    if (window() && (event_name.find("solve_") == 0 ||
                     event_name == "state_advance_before" ||
                     event_name == "state_advance_after" ||
                     event_name == "stepper_accept_after"))
    {
      const auto * solution = getNonlinearSystem(0).currentSolution();
      const auto field = event_name.find("solve_") == 0 ? "solve_return_vector" : event;
      if (solution)
        fingerprintVector(field, *solution);
      else
        _history_incomplete = true;
    }
    if (window())
      status();
  }
  catch (...)
  {
    _history_incomplete = true;
  }
}

void
CDPObservationProblem::computeResidual(const NumericVector<Number> & soln,
                                       NumericVector<Number> & residual,
                                       unsigned int sys_num)
{
  if (!_enabled)
  {
    FEProblem::computeResidual(soln, residual, sys_num);
    return;
  }
  CDPObservationFEContextGuard restore(_evaluation, _context, _global_iteration);
  CDPObservationContextGuard restore_point(_element, _qp, _call);
  ++_evaluation;
  _context = "residual";
  _global_iteration = -1;
  _element = 0;
  _qp = 0;
  _call = 0;
  if (window())
  {
    try
    {
      auto * solver = dynamic_cast<libMesh::PetscNonlinearSolver<Number> *>(
          getNonlinearSystem(sys_num).nonlinearSolver());
      PetscInt iteration = -1;
      Vec accepted = nullptr;
      if (!solver || SNESGetIterationNumber(solver->snes(), &iteration) ||
          SNESGetSolution(solver->snes(), &accepted))
        _history_incomplete = true;
      else
      {
        _global_iteration = iteration;
        const auto * vector = dynamic_cast<const libMesh::PetscVector<Number> *>(&soln);
        write("attempt-events",
              "residual_vector_identity",
              "identity",
              {},
              vector && vector->vec() == accepted ? "accepted_solution_vector" : "trial_vector");
      }
      fingerprintVector("residual_input_vector", soln);
      trialVector(soln);
    }
    catch (...)
    {
      _history_incomplete = true;
    }
  }
  try
  {
    FEProblem::computeResidual(soln, residual, sys_num);
  }
  catch (...)
  {
    write("attempt-events", "residual_exception", "completed", {}, "base residual threw");
    throw;
  }
  if (window())
    write("attempt-events", "residual_complete", "completed", {}, "no additional collective norm");
}

void
CDPObservationProblem::computeJacobian(const NumericVector<Number> & soln,
                                       SparseMatrix<Number> & jacobian,
                                       unsigned int sys_num)
{
  if (!_enabled)
  {
    FEProblem::computeJacobian(soln, jacobian, sys_num);
    return;
  }
  CDPObservationFEContextGuard restore(_evaluation, _context, _global_iteration);
  CDPObservationContextGuard restore_point(_element, _qp, _call);
  ++_evaluation;
  _context = "jacobian";
  _global_iteration = -1;
  _element = 0;
  _qp = 0;
  _call = 0;
  try
  {
    FEProblem::computeJacobian(soln, jacobian, sys_num);
  }
  catch (...)
  {
    write("attempt-events", "jacobian_exception", "completed", {}, "base jacobian threw");
    throw;
  }
  if (window())
    write("attempt-events", "jacobian_complete", "completed", {}, "base jacobian returned");
}

void
CDPObservationProblem::fingerprintVector(const char * field,
                                         const NumericVector<Number> & soln) noexcept
{
  try
  {
    const auto first = soln.first_local_index();
    const auto last = soln.last_local_index();
    CDPObservationFingerprint fingerprint(first, last);
    for (auto i = first; i < last; ++i)
      fingerprint.add(soln(i));
    std::ostringstream note;
    note << "fnv1a64:" << std::hex << std::setfill('0') << std::setw(16) << fingerprint.value()
         << "; rank-local candidate filter only, not proof of acceptance";
    CDPObservationContextGuard restore(_element, _qp, _call);
    _element = 0;
    _qp = 0;
    _call = 0;
    write("attempt-events",
          "vector_fingerprint",
          field,
          {static_cast<double>(first), static_cast<double>(last)},
          note.str());
  }
  catch (...)
  {
    _history_incomplete = true;
  }
}

void
CDPObservationProblem::trialVector(const NumericVector<Number> & soln)
{
  const auto & sys = getNonlinearSystem(0).system();
  const auto & vertices = getParam<std::vector<Real>>("observation_vertices");
  _element = 0;
  _qp = 0;
  for (const auto * node : mesh().getMesh().local_node_ptr_range())
  {
    bool selected = false;
    for (unsigned int j = 0; j < vertices.size(); j += 2)
      if (std::abs((*node)(0) - vertices[j]) < 1e-7 &&
          std::abs((*node)(1) - vertices[j + 1]) < 1e-7)
        selected = true;
    if (!selected)
      continue;
    for (unsigned int variable = 0; variable < sys.n_vars(); ++variable)
    {
      if (!node->n_dofs(sys.number(), variable))
        continue;
      const auto dof = node->dof_number(sys.number(), variable, 0);
      if (dof >= soln.first_local_index() && dof < soln.last_local_index())
        write("attempt-events",
              "residual_input",
              sys.variable_name(variable).c_str(),
              {static_cast<double>(node->id()),
               static_cast<double>(dof),
               soln(dof)},
              "node_dof_trial_value; vector identity is not alone proof of line-search acceptance");
    }
  }
}
