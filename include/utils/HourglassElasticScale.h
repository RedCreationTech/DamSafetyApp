#pragma once

namespace HourglassElasticScale
{
// Add only the elastic HHT history correction; damping is left in the original kernel.
template <typename Scalar>
Scalar correction(double scale, double alpha, double dt, const Scalar & current, double old)
{
  if (dt > 0)
    return (scale - 1.0) * ((1.0 + alpha) * current - alpha * old);
  return (scale - 1.0) * current;
}
}
