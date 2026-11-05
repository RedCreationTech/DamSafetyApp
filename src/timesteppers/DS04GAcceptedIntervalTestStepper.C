#include "DS04GAcceptedIntervalTestStepper.h"

#include "FEProblemBase.h"
#include "MooseMesh.h"

#include <iomanip>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
bool
readCsvRow(std::istream & input, std::vector<std::string> & row)
{
  row.clear();
  std::string field;
  bool quoted = false;
  bool read_any = false;
  char c;
  while (input.get(c))
  {
    read_any = true;
    if (quoted)
    {
      if (c == '"')
      {
        if (input.peek() == '"')
        {
          input.get(c);
          field += '"';
        }
        else
          quoted = false;
      }
      else
        field += c;
    }
    else if (c == '"' && field.empty())
      quoted = true;
    else if (c == ',')
    {
      row.push_back(field);
      field.clear();
    }
    else if (c == '\n')
    {
      row.push_back(field);
      return true;
    }
    else if (c != '\r')
      field += c;
  }
  if (!read_any || quoted)
    return false;
  row.push_back(field);
  return true;
}

std::vector<std::vector<std::string>>
readCsvTable(const std::string & path)
{
  std::ifstream input(path);
  if (!input)
    mooseError("Missing observer CSV: ", path);
  std::vector<std::vector<std::string>> rows;
  std::vector<std::string> row;
  while (readCsvRow(input, row))
  {
    if (row.size() != 26)
      mooseError("Observer CSV does not have 26 fields: ", path);
    rows.push_back(row);
  }
  if (!input.eof() || rows.empty())
    mooseError("Invalid or empty observer CSV: ", path);
  return rows;
}

std::vector<std::string>
numericValues(const std::vector<std::string> & row)
{
  std::stringstream stream(row[24]);
  std::vector<std::string> values;
  std::string value;
  while (stream >> value)
    values.push_back(value);
  return values;
}

std::string
materialCallKey(const std::vector<std::string> & row)
{
  return row[6] + ":" + row[8] + ":" + row[7] + ":" + row[10];
}
}

registerMooseObject("DamSafetyApp", DS04GAcceptedIntervalTestStepper);

InputParameters
DS04GAcceptedIntervalTestStepper::validParams()
{
  auto parameters = CDPObservationConstantDT::validParams();
  parameters.addClassDescription("TEST ONLY: one injected failure to verify real accepted-dt semantics");
  parameters.addParam<bool>(
      "check_observation", false, "Check the passive observer event stream after execution");
  parameters.addParam<bool>("check_material_observation",
                            false,
                            "Check real CDP, local/substep, swapBack and accepted/output records");
  return parameters;
}

DS04GAcceptedIntervalTestStepper::DS04GAcceptedIntervalTestStepper(const InputParameters & parameters)
  : CDPObservationConstantDT(parameters),
    _check_observation(getParam<bool>("check_observation")),
    _check_material_observation(getParam<bool>("check_material_observation")),
    _trace("accepted_interval_executioner.csv")
{
  _trace << std::setprecision(17) << "event,step,time,h,h_previous,last_accepted_h\n";
}

void
DS04GAcceptedIntervalTestStepper::record(const char * event)
{
  _trace << event << ',' << _t_step << ',' << _time << ',' << _dt << ','
         << _fe_problem.dtOld() << ',' << _last_accepted << '\n';
  if (!_trace)
    mooseError("Could not write accepted interval trace");
}

void
DS04GAcceptedIntervalTestStepper::step()
{
  if (comm().size() != 1 || _end_time > .01 || _accepted >= 8 ||
      _fe_problem.mesh().getMesh().n_elem() != 1)
    mooseError("Accepted interval tester is bounded to serial one element and at most .01s/8 steps");
  record("attempt");
  if (_accepted && _fe_problem.dtOld() != _last_accepted)
    mooseError("Previous accepted dt was changed by a rejected attempt");
  if (_cutback_accepted && _dt > _last_accepted)
    _after_cutback = true;
  CDPObservationConstantDT::step();
  if (_t_step == 2 && !_injected)
  {
    if (!_converged)
      mooseError("The bounded test solve unexpectedly failed before injection");
    _injected = true;
    _converged = false;
    ++_failure_count;
    record("injected_failure");
  }
}

void
DS04GAcceptedIntervalTestStepper::acceptStep()
{
  record("accept");
  if (_injected && _dt < _last_accepted)
    _cutback_accepted = true;
  _last_accepted = _dt;
  ++_accepted;
  CDPObservationConstantDT::acceptStep();
}

void
DS04GAcceptedIntervalTestStepper::rejectStep()
{
  const Real previous = _fe_problem.dtOld();
  CDPObservationConstantDT::rejectStep();
  if (_fe_problem.dtOld() != previous)
    mooseError("rejectStep modified accepted dt");
  record("reject_restore");
}

void
DS04GAcceptedIntervalTestStepper::verifyObservation() const
{
  std::ifstream input("accepted-output-observer-attempt-events-rank0.csv");
  if (!input)
    mooseError("Missing accepted-output observer event stream");
  std::vector<std::string> row;
  if (!readCsvRow(input, row))
    mooseError("Empty accepted-output observer event stream");
  const std::vector<std::string> expected_header = {
      "rank", "record", "time", "dt", "dt_old", "step", "attempt", "evaluation",
      "execution", "output", "call", "global_iteration", "context", "exec_flag", "element",
      "source_element", "qp", "partition", "N", "substep", "local_solve", "local_iteration",
      "event", "field", "values", "note"};
  if (row != expected_header)
    mooseError("Unexpected accepted-output observer schema");

  std::uint64_t previous_record = 0;
  std::uint64_t previous_attempt = 0;
  std::uint64_t previous_execution = 0;
  std::uint64_t previous_output = 0;
  std::uint64_t reject_attempt = 0;
  std::string reject_step;
  unsigned int attempt_enters = 0;
  unsigned int accepted_events = 0;
  unsigned int reject_before = 0;
  unsigned int reject_after = 0;
  unsigned int timestep_end_executions = 0;
  unsigned int timestep_end_outputs = 0;
  bool retry_found = false;
  while (readCsvRow(input, row))
  {
    if (row.size() != 26)
      mooseError("Observer event row does not have 26 fields");
    const auto record = std::stoull(row[1]);
    if (record <= previous_record || row[0] != "0")
      mooseError("Observer record/rank identity is not monotonic");
    previous_record = record;
    const auto & event = row[22];
    if (event == "attempt_enter")
    {
      ++attempt_enters;
      const auto attempt = std::stoull(row[6]);
      if (attempt != previous_attempt + 1)
        mooseError("Observer attempts do not include rejected retries");
      previous_attempt = attempt;
      if (reject_attempt && attempt == reject_attempt + 1 && row[5] == reject_step)
        retry_found = true;
    }
    else if (event == "accepted_end" || event == "stepper_accept_after")
    {
      if (row[12] != "accept" || row[7] != "0")
        mooseError("Accept event inherited an evaluation context");
      if (event == "accepted_end")
        ++accepted_events;
    }
    else if (event == "reject_before" || event == "reject_after")
    {
      if (row[12] != "reject" || row[7] != "0")
        mooseError("Reject event inherited an evaluation context");
      if (event == "reject_before")
      {
        ++reject_before;
        reject_attempt = std::stoull(row[6]);
        reject_step = row[5];
      }
      else
      {
        ++reject_after;
        if (std::stoull(row[6]) != reject_attempt)
          mooseError("Reject events do not share their attempt");
      }
    }
    else if (event == "execute_begin" || event == "execute_complete")
    {
      if (row[12] != "execute" || row[7] != "0")
        mooseError("Execute event inherited an evaluation context");
      const auto execution = std::stoull(row[8]);
      if (event == "execute_begin")
      {
        if (execution <= previous_execution)
          mooseError("Execution identity is not monotonic");
        previous_execution = execution;
        if (row[13] == "TIMESTEP_END")
          ++timestep_end_executions;
      }
      else if (execution != previous_execution)
        mooseError("Execute completion does not match its begin event");
    }
    else if (event == "output_begin" || event == "output_complete")
    {
      if (row[12] != "output" || row[7] != "0")
        mooseError("Output event inherited an evaluation context");
      const auto output = std::stoull(row[9]);
      if (event == "output_begin")
      {
        if (output <= previous_output)
          mooseError("Output identity is not monotonic");
        previous_output = output;
        if (row[13] == "TIMESTEP_END")
          ++timestep_end_outputs;
      }
      else if (output != previous_output)
        mooseError("Output completion does not match its begin event");
    }
  }
  if (!input.eof() || attempt_enters != _accepted + 1 || accepted_events != _accepted ||
      reject_before != 1 || reject_after != 1 || !retry_found ||
      timestep_end_executions < _accepted || timestep_end_outputs < _accepted)
    mooseError("Incomplete accepted/reject/execute/output observer sequence");
}

void
DS04GAcceptedIntervalTestStepper::verifyMaterialObservation() const
{
  const auto updates = readCsvTable("accepted-output-observer-qp-updates-rank0.csv");
  std::set<std::string> completed;
  std::set<std::string> enter_fields;
  std::set<std::string> success_fields;
  unsigned int enters = 0;
  unsigned int successes = 0;
  for (std::size_t i = 1; i < updates.size(); ++i)
  {
    const auto & row = updates[i];
    if (row[22] == "material_enter")
    {
      ++enters;
      enter_fields.insert(row[23]);
      const auto values = numericValues(row);
      if ((row[23] == "old16" && values.size() != 16) ||
          (row[23] == "rotation_row_major" && values.size() != 9) ||
          (row[23] == "dt_full_tangent" && values.size() != 2) ||
          ((row[23] == "old_elastic_strain" || row[23] == "raw_increment" ||
            row[23] == "old_total_strain" || row[23] == "new_total_strain") &&
           values.size() != 6))
        mooseError("Invalid material-entry observer payload");
    }
    else if (row[22] == "material_success")
    {
      ++successes;
      success_fields.insert(row[23]);
      completed.insert(materialCallKey(row));
      const auto values = numericValues(row);
      if ((row[23] == "new16" && values.size() != 16) ||
          (row[23] == "cauchy_stress" && values.size() != 6) ||
          (row[23] == "N_cutbacks_partitions_branch_plastic_final_local_residual" &&
           values.size() != 6))
        mooseError("Invalid material-exit observer payload");
    }
  }
  const std::set<std::string> expected_enter = {"old_elastic_strain",
                                                "raw_increment",
                                                "old_total_strain",
                                                "new_total_strain",
                                                "old16",
                                                "rotation_row_major",
                                                "dt_full_tangent"};
  const std::set<std::string> expected_success = {
      "cauchy_stress", "new16", "N_cutbacks_partitions_branch_plastic_final_local_residual"};
  if (!enters || enters != successes || enter_fields != expected_enter ||
      success_fields != expected_success)
    mooseError("Incomplete real-CDP material entry/exit records");

  const auto local = readCsvTable("accepted-output-observer-local-newton-rank0.csv");
  bool trial = false;
  for (std::size_t i = 1; i < local.size(); ++i)
    if (local[i][22] == "trial" && local[i][23] == "trial_stress" &&
        numericValues(local[i]).size() == 6 && std::stoull(local[i][20]) > 0)
      trial = true;
  if (!trial)
    mooseError("Missing real-CDP local trial record");

  const auto substeps = readCsvTable("accepted-output-observer-substep-events-rank0.csv");
  bool partitioned = false;
  for (std::size_t i = 1; i < substeps.size(); ++i)
    if (substeps[i][22] == "partition_succeeded" && std::stoull(substeps[i][18]) > 1)
      partitioned = true;
  if (!partitioned)
    mooseError("Missing real-CDP partitioned substep record");

  const auto swaps = readCsvTable("accepted-output-observer-swapback-events-rank0.csv");
  unsigned int swap_after = 0;
  unsigned int candidates = 0;
  for (std::size_t i = 1; i < swaps.size(); ++i)
  {
    if (swaps[i][22] == "swapback_after")
      ++swap_after;
    else if (swaps[i][22] == "swapback_candidate")
    {
      ++candidates;
      const auto values = numericValues(swaps[i]);
      if (values.size() != 4)
        mooseError("Invalid swapBack candidate key");
      const auto key = values[0] + ":" + values[1] + ":" + values[2] + ":" + values[3];
      if (!completed.count(key))
        mooseError("SwapBack candidate does not map to a completed material call");
    }
  }
  if (!swap_after || !candidates)
    mooseError("Missing swapBack state/candidate records");

  const auto history = readCsvTable("accepted-output-observer-history-events-rank0.csv");
  bool swap_history = false;
  for (std::size_t i = 1; i < history.size(); ++i)
    if (history[i][22] == "material_swapback_after")
      swap_history = true;
  if (!swap_history)
    mooseError("Missing post-swapBack material history");

  std::ifstream status_file("accepted-output-observer-status-rank0.json");
  std::stringstream status_buffer;
  status_buffer << status_file.rdbuf();
  const auto status = status_buffer.str();
  if (status.find("\"capture_truncated\":false") == std::string::npos ||
      status.find("\"io_failed\":false") == std::string::npos ||
      status.find("\"history_incomplete\":false") == std::string::npos ||
      status.find("\"swap_trace_incomplete\":false") == std::string::npos)
    mooseError("Real-CDP observer status is incomplete");
}

void
DS04GAcceptedIntervalTestStepper::postExecute()
{
  CDPObservationConstantDT::postExecute();
  if (!_injected || !_cutback_accepted || !_after_cutback || _accepted > 8 || !_accepted)
    mooseError("Incomplete or unbounded accepted-interval test");
  if (_check_observation)
    verifyObservation();
  if (_check_material_observation)
    verifyMaterialObservation();
  mooseInfo("DS04G accepted interval executioner check complete");
}
