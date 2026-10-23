#pragma once

#include "metaphysicl/raw_type.h"
#include <array>
#include <cmath>
#include <stdexcept>

namespace AqueductShell
{
template <typename T>
struct Stress
{
  std::array<T, 4> value; // xx, yy, zz, xy; transverse shear is section-elastic
  T plastic_strain;
  double zz_tangent;
};

// Three-dimensional J2 return with the source's two-row hardening table and
// constant extrapolation beyond its last plastic-strain value.
template <typename T>
Stress<T>
returnMap(const std::array<T, 4> & de,
          const std::array<double, 4> & old,
          double ep,
          double E,
          double nu,
          double yield,
          double ultimate,
          double limit)
{
  using std::sqrt;
  const double G = E / (2 * (1 + nu));
  const double K = E / (3 * (1 - 2 * nu));
  const double lambda = K - 2 * G / 3;
  const T trace = de[0] + de[1] + de[2];
  std::array<T, 4> trial;
  for (unsigned i = 0; i < 3; ++i)
    trial[i] = old[i] + 2 * G * de[i] + lambda * trace;
  trial[3] = old[3] + 2 * G * de[3];
  const T p = (trial[0] + trial[1] + trial[2]) / 3;
  const std::array<T, 4> s = {trial[0] - p, trial[1] - p, trial[2] - p, trial[3]};
  const T q2 = 1.5 * (s[0]*s[0] + s[1]*s[1] + s[2]*s[2] + 2*s[3]*s[3]);
  const double H = (ultimate - yield) / limit;
  const double strength = yield + H * std::min(ep, limit);
  if (MetaPhysicL::raw_value(q2) <= strength * strength)
    return {trial, T(ep), lambda + 2 * G};
  const T q = sqrt(q2);
  T dl = (q - strength) / (3 * G + (ep < limit ? H : 0));
  double slope = ep < limit ? H : 0;
  if (ep + MetaPhysicL::raw_value(dl) > limit)
  {
    dl = (q - ultimate) / (3 * G);
    slope = 0;
  }
  const T factor = 1 - 3 * G * dl / q;
  Stress<T> out;
  for (unsigned i = 0; i < 3; ++i)
    out.value[i] = p + factor * s[i];
  out.value[3] = factor * s[3];
  out.plastic_strain = ep + dl;
  const double qr = MetaPhysicL::raw_value(q);
  const double sz = MetaPhysicL::raw_value(s[2]);
  const double qz = 3 * G * sz / qr;
  const double df = -3 * G * (qz / ((3 * G + slope) * qr) -
                              MetaPhysicL::raw_value(dl) * qz / (qr * qr));
  out.zz_tangent = K + MetaPhysicL::raw_value(factor) * 4 * G / 3 + sz * df;
  return out;
}

template <typename T>
Stress<T>
planeStress(const std::array<T, 3> & increment, // tensor strain xx, yy, xy
            const std::array<double, 4> & old,
            double ep,
            double E,
            double nu,
            double yield,
            double ultimate,
            double limit)
{
  if (!(E > 0 && nu > -1 && nu < .5 && yield > 0 && ultimate >= yield && limit > 0))
    throw std::invalid_argument("invalid aqueduct shell material");
  T ezz = -nu / (1 - nu) * (increment[0] + increment[1]);
  for (unsigned it = 0; it < 40; ++it)
  {
    auto state = returnMap<T>({increment[0], increment[1], ezz, increment[2]},
                              old, ep, E, nu, yield, ultimate, limit);
    const bool converged = std::abs(MetaPhysicL::raw_value(state.value[2])) < yield * 1e-11;
    // Correct even a zero-valued residual once so AD also satisfies d(sigma33)=0.
    ezz -= state.value[2] / state.zz_tangent;
    if (converged)
      return returnMap<T>({increment[0], increment[1], ezz, increment[2]},
                          old, ep, E, nu, yield, ultimate, limit);
  }
  throw std::runtime_error("aqueduct shell plane-stress return did not converge");
}

struct Moments
{
  double zero, first, second;
};

// Offset is the reference-surface position divided by thickness: SPOS = +0.5.
inline Moments moments(double thickness, double offset)
{
  if (!(thickness > 0) || !std::isfinite(offset))
    throw std::invalid_argument("invalid aqueduct shell section");
  return {thickness, -offset * thickness * thickness,
          thickness * thickness * thickness * (1. / 12 + offset * offset)};
}

inline std::array<double, 3> layerPositions(double thickness, double offset)
{
  return {thickness * (-.5 - offset), -thickness * offset,
          thickness * (.5 - offset)};
}

constexpr std::array<double, 3> layerWeights = {1. / 6, 4. / 6, 1. / 6};
}
