#pragma once
#include <cmath>
#include <stdexcept>

namespace AqueductHHTScheme
{
inline double massScale(double alpha)
{
  if (!std::isfinite(alpha) || alpha < -.5 || alpha > 0.)
    throw std::invalid_argument("HHT alpha must be in [-.5,0]");
  return 1. / (1. + alpha);
}
inline double beta(double alpha) { massScale(alpha); return .25 * (1.-alpha) * (1.-alpha); }
inline double gamma(double alpha) { massScale(alpha); return .5-alpha; }
template <typename T>
T normalizedResidual(const T & inertia, const T & force, const T & old_force, double alpha)
{
  return massScale(alpha) * inertia + force - alpha * massScale(alpha) * old_force;
}
}
