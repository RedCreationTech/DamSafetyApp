#include "AbsoluteGroundMotionAux.h"
#include "Function.h"
registerMooseObject("DamSafetyApp", AbsoluteGroundMotionAux);
InputParameters AbsoluteGroundMotionAux::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("Explicit absolute output: relative solution plus the same ground function used in inertia.");
  params.addRequiredCoupledVar("relative", "Relative displacement, velocity or acceleration");
  params.addRequiredParam<FunctionName>("ground", "Matching analytic ground component");
  return params;
}
AbsoluteGroundMotionAux::AbsoluteGroundMotionAux(const InputParameters & parameters)
  : AuxKernel(parameters), _relative(coupledValue("relative")), _ground(getFunction("ground"))
{
  if (!isNodal()) mooseError(name(), " requires a nodal auxiliary variable");
}
Real AbsoluteGroundMotionAux::computeValue()
{ return _relative[_qp] + _ground.value(_t, *_current_node); }
