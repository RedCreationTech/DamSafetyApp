#pragma once
#include "RankFourTensor.h"
#include "RankTwoTensor.h"
#include <stdexcept>

namespace DS04GFalseRotation
{
inline RankTwoTensor increment(const RankTwoTensor & Fhat)
{
  const auto G = Fhat.inverse();
  const auto C = G * G.transpose() - RankTwoTensor(RankTwoTensor::initIdentity);
  return -.5 * C + .25 * C * C;
}

inline RankTwoTensor incrementDerivative(const RankTwoTensor & Fhat, const RankTwoTensor & dFhat)
{
  const auto G = Fhat.inverse();
  const auto C = G * G.transpose() - RankTwoTensor(RankTwoTensor::initIdentity);
  const auto dG = -G * dFhat * G;
  const auto dC = dG * G.transpose() + G * dG.transpose();
  return -.5 * dC + .25 * (dC * C + C * dC);
}

struct Parts
{
  Real material;
  Real geometry;
};

inline Parts contract(const RankTwoTensor & F, const RankTwoTensor & Fold,
                      const RankFourTensor & D, const RankTwoTensor & weighted_stress,
                      Real current_weight, unsigned component, unsigned direction,
                      const RealGradient & test, const RealGradient & trial)
{
  if (F.det() <= 0 || Fold.det() <= 0)
    throw std::invalid_argument("SP1 false-rotation tangent requires positive F/Fold determinants");
  const auto Fhat = F * Fold.inverse();
  RankTwoTensor L;
  for (unsigned j = 0; j < 3; ++j)
    L(direction, j) = trial(j);
  const auto dT = incrementDerivative(Fhat, L * Fhat);
  return {current_weight * (D * dT).row(component) * test,
          weighted_stress.row(component) *
              (test * trial(direction) - trial * test(direction))};
}
}
