#pragma once

#include "MooseTypes.h"
#include "libmesh/point.h"
#include <array>

namespace UniformHex8
{
struct Operators
{
  Real volume = 0.0;
  std::array<RealVectorValue, 8> gradient;
  std::array<std::array<Real, 4>, 8> hourglass;
  Real gradient_norm_squared = 0.0;
  Real hourglassStiffness(unsigned int i, unsigned int j, Real shear_modulus,
                        Real factor = 0.005) const;
};

// Average gradients are integrated on the reference geometry, independently
// of the single constitutive integration point. No material is evaluated here.
Operators compute(const std::array<libMesh::Point, 8> & points);
}
