#pragma once
#include <stdexcept>

namespace StressHistoryInterval
{
// Multiplies (stress_older - stress_old). The current-stress tangent is unchanged.
inline double correction(double alpha, double zeta, double dt, double previous_dt, int step)
{
  if (dt <= 0 || step <= 1 || alpha == 0 || zeta == 0)
    return 0;
  if (previous_dt <= 0)
    throw std::invalid_argument("Previous accepted stress interval must be positive");
  return alpha * zeta * (1 / previous_dt - 1 / dt);
}
}
