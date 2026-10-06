#include "UniformHex8.h"
#include "RankTwoTensor.h"
#include <cmath>
#include <stdexcept>

namespace UniformHex8
{
Operators compute(const std::array<libMesh::Point, 8> & points)
{
  constexpr int sign[8][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
                             {-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
  Operators result;
  libMesh::Point center;
  for (const auto & point : points)
    center += point / 8.0;
  const Real gauss = 1.0 / std::sqrt(3.0);
  for (int gx : {-1, 1})
    for (int gy : {-1, 1})
      for (int gz : {-1, 1})
      {
        const Real xi[3] = {gx * gauss, gy * gauss, gz * gauss};
        std::array<RealVectorValue, 8> derivative;
        RankTwoTensor jacobian;
        for (unsigned int i = 0; i < 8; ++i)
          for (unsigned int a = 0; a < 3; ++a)
          {
            Real value = 0.125 * sign[i][a];
            for (unsigned int b = 0; b < 3; ++b)
              if (a != b)
                value *= 1.0 + sign[i][b] * xi[b];
            derivative[i](a) = value;
            for (unsigned int d = 0; d < 3; ++d)
              jacobian(d, a) += (points[i](d) - center(d)) * value;
          }
        const Real determinant = jacobian.det();
        if (!std::isfinite(determinant) || determinant <= 0.0)
          throw std::invalid_argument("UniformHex8: non-positive reference Jacobian");
        const auto inverse = jacobian.inverse();
        result.volume += determinant;
        for (unsigned int i = 0; i < 8; ++i)
          for (unsigned int d = 0; d < 3; ++d)
            for (unsigned int a = 0; a < 3; ++a)
              result.gradient[i](d) += determinant * inverse(a, d) * derivative[i](a);
      }
  for (auto & gradient : result.gradient)
  {
    gradient /= result.volume;
    result.gradient_norm_squared += gradient.norm_sq();
  }
  // Flanagan-Belytschko modes corrected to annihilate every affine physical
  // displacement, including rigid translations/rotations on distorted bricks.
  for (unsigned int mode = 0; mode < 4; ++mode)
  {
    Real base[8];
    RealVectorValue moment;
    for (unsigned int j = 0; j < 8; ++j)
    {
      base[j] = mode == 0 ? sign[j][1] * sign[j][2] :
                mode == 1 ? sign[j][0] * sign[j][2] :
                mode == 2 ? sign[j][0] * sign[j][1] :
                            sign[j][0] * sign[j][1] * sign[j][2];
      moment += (points[j] - center) * base[j];
    }
    for (unsigned int i = 0; i < 8; ++i)
      result.hourglass[i][mode] = base[i] - result.gradient[i] * moment;
  }
  return result;
}

Real Operators::hourglassStiffness(unsigned int i, unsigned int j, Real shear_modulus,
                                 Real factor) const
{
  if (i >= 8 || j >= 8 || shear_modulus <= 0 || factor < 0)
    throw std::invalid_argument("UniformHex8: invalid hourglass stiffness input");
  Real product = 0.0;
  for (unsigned int mode = 0; mode < 4; ++mode)
    product += hourglass[i][mode] * hourglass[j][mode];
  return factor * shear_modulus * volume * gradient_norm_squared * product;
}
}
