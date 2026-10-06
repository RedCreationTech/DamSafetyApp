#include "SmallSlidingContactTraction.h"
registerMooseObject("DamSafetyApp", SmallSlidingContactTraction);
InputParameters SmallSlidingContactTraction::validParams()
{
  auto p=ADMortarLagrangeConstraint::validParams();
  p.addRequiredRangeCheckedParam<unsigned int>("component","component < 3","Force component.");
  p.addRequiredParam<UserObjectName>("contact","Small sliding cohesive contact object.");
  p.set<bool>("compute_lm_residuals")=false;
  p.set<bool>("use_displaced_mesh")=false;
  p.set<bool>("interpolate_normals")=false;
  return p;
}
SmallSlidingContactTraction::SmallSlidingContactTraction(const InputParameters & p)
 : ADMortarLagrangeConstraint(p),_component(getParam<unsigned int>("component")),
   _contact(getUserObject<SmallSlidingCohesiveContact>("contact")) {}
ADReal SmallSlidingContactTraction::computeQpResidual(Moose::MortarType type)
{
  if(type==Moose::MortarType::Secondary) return _test_secondary[_i][_qp]*_contact.traction(_component)[_qp];
  if(type==Moose::MortarType::Primary) return -_test_primary[_i][_qp]*_contact.traction(_component)[_qp];
  return 0;
}
