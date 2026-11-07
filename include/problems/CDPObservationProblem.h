#pragma once

#include "FEProblem.h"
#include "CDPObservation.h"
#include "CDPObservationSwapTracker.h"

#include <fstream>
#include <map>
#include <petscsnes.h>
#include <set>
#include <string>

class CDPObservationProblem : public FEProblem, public CDPObservation::Sink
{
public:
  static InputParameters validParams();
  CDPObservationProblem(const InputParameters & parameters);
  ~CDPObservationProblem() override;

  void initialSetup() override;
  void advanceState() override;
  void swapBackMaterials(const THREAD_ID tid) override;
  void execute(const ExecFlagType & exec_type) override;
  void outputStep(ExecFlagType type) override;
  void computeResidual(const NumericVector<Number> & soln,
                       NumericVector<Number> & residual,
                       unsigned int sys_num) override;
  void computeJacobian(const NumericVector<Number> & soln,
                       SparseMatrix<Number> & jacobian,
                       unsigned int sys_num) override;

  bool enterPoint(std::uint64_t element, unsigned int qp, bool new_call) noexcept override;
  void materialCompleted() noexcept override;
  void failed() noexcept override { _io_failed = true; }
  void write(const char * stream,
             const char * event,
             const char * name,
             const std::vector<double> & values,
             const std::string & note) noexcept override;

  bool observationEnabled() const { return _enabled; }
  void attemptEvent(const char * event) noexcept;

private:
  bool window() const;
  void history(const char * event);
  void trialVector(const NumericVector<Number> & soln);
  void fingerprintVector(const char * field, const NumericVector<Number> & soln) noexcept;
  void setExecFlag(const ExecFlagType & type) noexcept;
  void status() noexcept;
  void residualGroups(const char * event) noexcept;
  std::vector<double> acceptedFingerprint();
  void restoreKernelEnable(const std::map<std::string, bool> & enabled) noexcept;

  const bool _enabled;
  const bool _trace_swap;
  const Real _start;
  const Real _end;
  const std::string _prefix;
  const std::uint64_t _limit;

  CDPObservationSwapTracker _swap_tracker{64};
  std::uint64_t _swap_events = 0;
  std::uint64_t _swap_empty = 0;
  bool _swap_incomplete = false;

  std::map<std::uint64_t, unsigned int> _targets;
  std::map<std::string, std::ofstream> _streams;
  std::map<std::string, std::uint64_t> _counts;
  double _write_microseconds = 0.0;

  std::map<std::pair<std::uint64_t, unsigned int>, CDPObservationSwapTracker::Call> _last_call;
  std::uint64_t _next_call = 0;
  std::uint64_t _attempt = 0;
  std::uint64_t _evaluation = 0;
  std::uint64_t _execution = 0;
  std::uint64_t _output = 0;
  std::uint64_t _call = 0;
  std::uint64_t _record = 0;
  std::uint64_t _bytes = 0;
  std::uint64_t _element = 0;
  unsigned int _qp = 0;
  const char * _context = "event";
  int _global_iteration = -1;
  std::string _exec_flag = "NONE";
  bool _truncated = false;
  bool _io_failed = false;
  bool _history_incomplete = false;

  const bool _res_groups;
  const std::map<std::string, std::vector<std::string>> _res_group_terms;
  const std::vector<unsigned int> _res_extra_nodes;
  const BoundaryName _res_boundary;
  std::vector<dof_id_type> _res_target_nodes;
  std::set<dof_id_type> _res_base_nodes;
  std::uint64_t _res_rows = 0;
  bool _res_contaminated = false;
};
