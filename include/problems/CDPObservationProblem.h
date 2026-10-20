#pragma once
#include "FEProblem.h"
#include "CDPObservation.h"
#include <fstream>
#include <map>
#include <petscsnes.h>

/** Passive DS04G recorder. No additional residual, tangent, or material evaluations. */
class CDPObservationProblem : public FEProblem, public CDPObservation::Sink
{
public:
  static InputParameters validParams();
  CDPObservationProblem(const InputParameters &);
  ~CDPObservationProblem() override;
  void initialSetup() override;
  void computeResidual(const NumericVector<Number> &, NumericVector<Number> &, unsigned int) override;
  void computeJacobian(const NumericVector<Number> &, SparseMatrix<Number> &, unsigned int) override;
  bool enterPoint(std::uint64_t, unsigned int, bool) override;
  void write(const char *, const char *, const char *, const std::vector<double> &,
             const std::string &) noexcept override;
  void attemptEvent(const char *);
  static PetscErrorCode monitor(SNES, PetscInt, PetscReal, void *);
private:
  bool window() const;
  void history(const char *);
  void status() noexcept;
  const bool _enabled;
  const Real _start, _end, _offset;
  const std::string _prefix;
  const std::uint64_t _limit;
  std::map<std::uint64_t, unsigned int> _targets;
  std::map<std::string, std::ofstream> _streams;
  std::map<std::string, std::uint64_t> _counts;
  double _write_microseconds = 0.;
  std::map<std::pair<std::uint64_t, unsigned int>, std::pair<std::uint64_t, std::uint64_t>> _last_call;
  std::uint64_t _next_call = 0;
  std::uint64_t _attempt = 0, _evaluation = 0, _call = 0, _record = 0, _bytes = 0;
  std::uint64_t _element = 0;
  unsigned int _qp = 0;
  std::string _kind = "event";
  bool _truncated = false, _io_failed = false, _history_incomplete = false;
};
