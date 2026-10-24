#pragma once
#include <array>
#include <cmath>

namespace AqueductFriction
{
template <typename T> struct State
{
  std::array<T, 2> elastic;
  std::array<T, 2> traction;
};

// Constant-mu Abaqus elastic predictor and circular return, with pressure-dependent stick stiffness.
template <typename T>
State<T> update(const std::array<double, 2> & old,
                const std::array<T, 2> & increment,
                const T & pressure,
                double mu,
                double critical_slip)
{
  State<T> out{};
  if (pressure <= 0.)
    return out;
  out.elastic = {old[0] + increment[0], old[1] + increment[1]};
  const T norm2 = out.elastic[0] * out.elastic[0] + out.elastic[1] * out.elastic[1];
  if (norm2 > critical_slip * critical_slip)
  {
    using std::sqrt;
    const T scale = critical_slip / sqrt(norm2);
    for (auto & e : out.elastic)
      e *= scale;
  }
  for (unsigned j = 0; j < 2; ++j)
    out.traction[j] = mu * pressure * out.elastic[j] / critical_slip;
  return out;
}
}
