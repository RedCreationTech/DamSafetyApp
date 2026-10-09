#pragma once
#include <stdexcept>

namespace ElasticVelocityRate
{
// Accepted u/v/a are immutable during a trial, including retries with a different dt.
template <typename T>
T trialVelocity(const T & u, const T & old_u, const T & old_v, const T & old_a,
                double dt, double beta, double gamma)
{
  if (!(dt > 0) || !(beta > 0))
    throw std::invalid_argument("ElasticVelocityRate requires positive dt and beta");
  const T a = ((u - old_u) / (dt * dt) - old_v / dt - old_a * (0.5 - beta)) / beta;
  return old_v + old_a * (dt * (1 - gamma)) + a * (gamma * dt);
}
inline double displacementDerivative(double dt, double beta, double gamma)
{
  if (!(dt > 0) || !(beta > 0))
    throw std::invalid_argument("ElasticVelocityRate requires positive dt and beta");
  return gamma / (beta * dt);
}
}
