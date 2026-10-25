#pragma once

#include "AbaqusCDPSubstepIntegrator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <limits>
#include <tuple>
#include <utility>

/** One successful complete result per material point, scoped to one assembly epoch.
 * Each material/thread owns its cache and immutable integrator configuration.
 * Hashes and numerical tolerances are deliberately not used for equality.
 */
class AbaqusCDPExactResultCache
{
public:
  using Result = AbaqusCDPSubstepIntegrator::LinearizedResult;
  struct Point
  {
    std::uint64_t first;
    unsigned int second, side;
    bool operator==(const Point & p) const
    { return std::tie(first,second,side)==std::tie(p.first,p.second,p.side); }
    bool operator<(const Point & p) const
    { return std::tie(first,second,side)<std::tie(p.first,p.second,p.side); }
  };
  struct Key
  {
    Point point;
    std::array<double, 29> input;
    bool finite() const
    {
      for (double value : input)
        if (!std::isfinite(value))
          return false;
      return true;
    }
    bool operator==(const Key & other) const
    {
      return point == other.point && finite() && other.finite() &&
             std::memcmp(input.data(), other.input.data(), sizeof(double) * input.size()) == 0;
    }
  };
  struct Counters
  {
    std::uint64_t requests = 0, potential_hits = 0, reused_hits = 0;
    std::uint64_t stores = 0, capacity_skips = 0, invalidations = 0;
    std::size_t peak_entries = 0;
  };

  explicit AbaqusCDPExactResultCache(std::size_t capacity) : _capacity(capacity) {}

  static Key key(std::uint64_t element,
                 unsigned int qp,
                 const AbaqusCDPSubstepIntegrator::SymmetricTensor & old_strain,
                 const AbaqusCDPSubstepIntegrator::SymmetricTensor & new_strain,
                 double dt,
                 const AbaqusCDPSubstepIntegrator::State & old,
                 unsigned int side=std::numeric_limits<unsigned int>::max())
  {
    Key result{{element, qp, side}, {}};
    std::size_t k = 0;
    for (double v : old_strain) result.input[k++] = v;
    for (double v : new_strain) result.input[k++] = v;
    result.input[k++] = dt;
    for (double v : old.backbone.plastic_strain) result.input[k++] = v;
    result.input[k++] = old.backbone.tensile_equivalent_plastic_strain;
    result.input[k++] = old.backbone.compressive_equivalent_plastic_strain;
    for (double v : old.viscous_plastic_strain) result.input[k++] = v;
    result.input[k++] = old.viscous_tension_damage;
    result.input[k++] = old.viscous_compression_damage;
    return result;
  }

  const Result * find(const Key & key, bool reuse)
  {
    ++_counters.requests;
    const auto it = _entries.find(key.point);
    if (it == _entries.end() || !(it->second.key == key))
      return nullptr;
    ++_counters.potential_hits;
    if (reuse)
      ++_counters.reused_hits;
    return &it->second.result;
  }

  void store(const Key & key, const Result & result)
  {
    if (!key.finite())
      return;
    auto it = _entries.find(key.point);
    if (it != _entries.end())
      it->second = {key, result};
    else if (_entries.size() < _capacity)
      _entries.emplace(key.point, Entry{key, result});
    else
    {
      ++_counters.capacity_skips;
      return;
    }
    ++_counters.stores;
    _counters.peak_entries = std::max(_counters.peak_entries, _entries.size());
  }

  void invalidate()
  {
    _entries.clear();
    ++_counters.invalidations;
  }
  const Counters & counters() const { return _counters; }
  std::size_t entries() const { return _entries.size(); }
  std::size_t entryPayloadBytes() const { return sizeof(Entry); }

private:
  struct Entry { Key key; Result result; };
  const std::size_t _capacity;
  std::map<Point, Entry> _entries;
  Counters _counters;
};
