#pragma once
#include <cstdint>

/** Restore recorder-only context even if an observer read/allocation fails. */
class CDPObservationContextGuard
{
public:
  CDPObservationContextGuard(std::uint64_t & element, unsigned int & qp, std::uint64_t & call)
    : _element(element), _qp(qp), _call(call), _old_element(element), _old_qp(qp), _old_call(call) {}
  ~CDPObservationContextGuard() noexcept { _element=_old_element; _qp=_old_qp; _call=_old_call; }
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
