#include "SmallSlidingContactAux.h"
registerMooseObject("DamSafetyApp",SmallSlidingContactAux);
InputParameters SmallSlidingContactAux::validParams()
{
  auto p=AuxKernel::validParams();
  p.addRequiredParam<UserObjectName>("contact","Small-sliding contact object.");
  p.addRequiredParam<MooseEnum>("quantity",MooseEnum("damage=0 gap=1 pressure=2 jump_n=3 jump_s=4 jump_t=5 traction_n=6 traction_s=7 traction_t=8"),"Contact quantity; tractions are residual-sign local stresses.");
  p.set<bool>("use_displaced_mesh")=false;
  return p;
}
SmallSlidingContactAux::SmallSlidingContactAux(const InputParameters & p)
 : AuxKernel(p),_contact(getUserObject<SmallSlidingCohesiveContact>("contact")),_quantity(getParam<MooseEnum>("quantity"))
{
  if(!isNodal())mooseError("Contact output must use a nodal auxiliary field.");
}
Real SmallSlidingContactAux::computeValue(){return _contact.nodeValue(_current_node,_quantity);}
