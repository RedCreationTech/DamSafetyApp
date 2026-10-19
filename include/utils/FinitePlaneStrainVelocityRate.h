#pragma once

#include "ElasticVelocityRate.h"
#include "RankTwoTensor.h"
#include <cmath>
#include <stdexcept>

namespace FinitePlaneStrainVelocityRate
{
template <typename T>
RankTwoTensorTempl<T> deformation(const RankTwoTensorTempl<T> & gradient)
{
  auto f = gradient;
  f.addIa(1.0);
  if (!(MetaPhysicL::raw_value(f.det()) > 0.0))
    throw std::invalid_argument("Finite Newmark damping requires positive det(F)");
  return f;
}

template <typename T>
RankTwoTensorTempl<T> stress(const RankTwoTensorTempl<T> & velocity_gradient,
                           const RankTwoTensorTempl<T> & f,
                           Real youngs, Real poisson)
{
  const auto l = velocity_gradient * f.inverse();
  auto d = (l + l.transpose()) * 0.5;
  auto result = d * (youngs / (1.0 + poisson));
  result.addIa(d.trace() * (youngs * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson))));
  return result;
}

template <typename T>
RankTwoTensorTempl<T> nominalStress(
    const RankTwoTensorTempl<T> & gradient,
    const RankTwoTensor & old_gradient,
    const RankTwoTensor & old_velocity_gradient,
    const RankTwoTensor & old_acceleration_gradient,
    Real dt, Real beta, Real gamma, Real alpha, Real zeta, Real youngs, Real poisson)
{
  if (zeta == 0.0 || dt <= 0.0)
    return RankTwoTensorTempl<T>();
  const auto f = deformation(gradient);
  const auto old_f = deformation(old_gradient);
  RankTwoTensorTempl<T> velocity_gradient;
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      velocity_gradient(i,j) = ElasticVelocityRate::trialVelocity(
          gradient(i,j), T(old_gradient(i,j)), T(old_velocity_gradient(i,j)),
          T(old_acceleration_gradient(i,j)), dt, beta, gamma);
  const auto current = stress(velocity_gradient, f, youngs, poisson);
  const auto old = stress(old_velocity_gradient, old_f, youngs, poisson);
  // Keep upstream HHT Cauchy-history weighting, then pull back the current weak form.
  const RankTwoTensorTempl<T> weighted = (current * (1.0 + alpha) - old * alpha) * zeta;
  return (weighted * f.inverse().transpose()) * f.det();
}
}
