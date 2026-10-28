#pragma once
#include "Gravity.h"
#include <array>

/** SP1-only optional geometric Jacobian; the residual is inherited unchanged. */
class DS04GSP1Gravity : public Gravity
{
public:
  static InputParameters validParams();
  DS04GSP1Gravity(const InputParameters &);

protected:
  void initialSetup() override;
  Real computeQpJacobian() override;
  Real computeQpOffDiagJacobian(unsigned) override;
  Real tangent(unsigned direction);
  const std::array<unsigned, 2> _displacement_vars;
  bool _complete = false;
};
