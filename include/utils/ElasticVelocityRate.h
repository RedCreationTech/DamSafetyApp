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

// Only for an independently audited observable plane-stress rate space.
// This eliminates rate_zz, never the original static thickness unknown.
template <typename Fourth, typename Second>
Second closePlaneStressRate(const Fourth & elasticity, Second rate)
{
  if (!(elasticity(2,2,2,2) > 0))
    throw std::invalid_argument("Elastic plane-stress rate requires positive Czzzz");
  rate(2,2) = 0;
  rate(2,2) = -(elasticity * rate)(2,2) / elasticity(2,2,2,2);
  return rate;
}

template <typename Fourth>
Fourth planeStressRateTangent(const Fourth & elasticity)
{
  if (!(elasticity(2,2,2,2) > 0))
    throw std::invalid_argument("Elastic plane-stress rate requires positive Czzzz");
  Fourth tangent = elasticity;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int l = 0; l < 3; ++l)
          tangent(i,j,k,l) -= elasticity(i,j,2,2) * elasticity(2,2,k,l) /
                              elasticity(2,2,2,2);
  return tangent;
}
}
