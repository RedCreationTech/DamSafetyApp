#pragma once
#include <cmath>
#include <stdexcept>
#include <limits>

namespace AqueductOutputSchedule
{
// Mark-anchored output never introduces a solver synchronization point.
inline long long slot(double time,double origin,double end,double interval,double tolerance)
{
  if(!std::isfinite(time) || !std::isfinite(origin) || !std::isfinite(end) ||
     !std::isfinite(interval) || !std::isfinite(tolerance) || end<origin ||
     interval<=0. || tolerance<0. || tolerance>=interval*.5)
    throw std::invalid_argument("Invalid aqueduct output schedule");
  if(time<origin-tolerance || time>end+tolerance)return -1;
  const double mark=std::floor((time-origin+tolerance)/interval);
  if(mark>static_cast<double>(std::numeric_limits<long long>::max()/2))
    throw std::invalid_argument("Aqueduct output schedule index overflow");
  return static_cast<long long>(mark);
}
}
