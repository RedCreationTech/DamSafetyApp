#pragma once

#include <cstdint>

class CDPObservationContextGuard
{
public:
  CDPObservationContextGuard(std::uint64_t & element, unsigned int & qp, std::uint64_t & call)
    : _element(element),
      _qp(qp),
      _call(call),
      _old_element(element),
      _old_qp(qp),
      _old_call(call)
  {
  }

  ~CDPObservationContextGuard() noexcept
  {
    _element = _old_element;
    _qp = _old_qp;
    _call = _old_call;
  }

  CDPObservationContextGuard(const CDPObservationContextGuard &) = delete;
  CDPObservationContextGuard & operator=(const CDPObservationContextGuard &) = delete;

private:
  std::uint64_t & _element;
  unsigned int & _qp;
  std::uint64_t & _call;
  const std::uint64_t _old_element;
  const unsigned int _old_qp;
  const std::uint64_t _old_call;
};

class CDPObservationFEContextGuard
{
public:
  CDPObservationFEContextGuard(std::uint64_t & evaluation,
                               const char *& context,
                               int & global_iteration) noexcept
    : _evaluation(evaluation),
      _context(context),
      _global_iteration(global_iteration),
      _old_evaluation(evaluation),
      _old_context(context),
      _old_global_iteration(global_iteration)
  {
  }

  ~CDPObservationFEContextGuard() noexcept
  {
    _evaluation = _old_evaluation;
    _context = _old_context;
    _global_iteration = _old_global_iteration;
  }

  CDPObservationFEContextGuard(const CDPObservationFEContextGuard &) = delete;
  CDPObservationFEContextGuard & operator=(const CDPObservationFEContextGuard &) = delete;

private:
  std::uint64_t & _evaluation;
  const char *& _context;
  int & _global_iteration;
  const std::uint64_t _old_evaluation;
  const char * const _old_context;
  const int _old_global_iteration;
};
