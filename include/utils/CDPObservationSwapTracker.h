#pragma once

#include <cstdint>
#include <map>
#include <utility>
#include <vector>

class CDPObservationSwapTracker
{
public:
  struct Call
  {
    std::uint64_t attempt;
    std::uint64_t execution;
    std::uint64_t evaluation;
    std::uint64_t call;
  };

  struct Batch
  {
    std::vector<Call> calls;
    bool overflow = false;
  };

  explicit CDPObservationSwapTracker(unsigned int limit) : _limit(limit) {}

  void completed(std::uint64_t element,
                 std::uint64_t attempt,
                 std::uint64_t execution,
                 std::uint64_t evaluation,
                 std::uint64_t call)
  {
    auto & batch = _pending[element];
    if (batch.calls.size() < _limit)
      batch.calls.push_back({attempt, execution, evaluation, call});
    else
      batch.overflow = true;
  }

  Batch take(std::uint64_t element)
  {
    auto it = _pending.find(element);
    if (it == _pending.end())
      return {};
    auto out = std::move(it->second);
    _pending.erase(it);
    return out;
  }

private:
  unsigned int _limit;
  std::map<std::uint64_t, Batch> _pending;
};
