#pragma once

// Only the missing XY shear terms; the upstream normal terms and residual stay intact.
namespace PlaneStressCoupling
{
template <typename Tensor, typename Gradient, typename Scalar>
auto displacementShearColumn(const Tensor & C, unsigned int component,
                             const Gradient & gradient_test, const Scalar & phi, double alpha)
{
  return (1. + alpha) * C(component, 1 - component, 2, 2) *
         gradient_test(1 - component) * phi;
}

template <typename Tensor, typename Scalar, typename Gradient>
auto thicknessShearRow(const Tensor & C, unsigned int component,
                       const Scalar & test, const Gradient & gradient_phi)
{
  return .5 * (C(2, 2, component, 1 - component) + C(2, 2, 1 - component, component)) *
         test * gradient_phi(1 - component);
}
}
