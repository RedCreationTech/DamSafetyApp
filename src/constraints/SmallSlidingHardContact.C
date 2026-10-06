#include "SmallSlidingHardContact.h"
#include "MooseVariable.h"
#include "SystemBase.h"
#include "Assembly.h"
registerMooseObject("DamSafetyApp", SmallSlidingHardContact);
InputParameters SmallSlidingHardContact::validParams()
{
  auto p=ADMortarConstraint::validParams();
  p.addClassDescription("Nodal mortar HARD contact complementarity on fixed reference geometry.");
  p.addRequiredParam<UserObjectName>("contact","Small sliding cohesive contact object.");
  p.addRangeCheckedParam<Real>("c",1e12,"c > 0","Pressure/gap balance, Pa/m; not penalty compliance.");
  p.set<bool>("use_displaced_mesh")=false;
  p.set<bool>("interpolate_normals")=false;
  p.set<bool>("compute_primal_residuals")=false;
  return p;
}
SmallSlidingHardContact::SmallSlidingHardContact(const InputParameters & p)
 : ADMortarConstraint(p),_c(getParam<Real>("c")),
   _contact(getUserObject<SmallSlidingCohesiveContact>("contact")) {}
void SmallSlidingHardContact::post()
{
  for(const auto & entry:_contact.weightedGaps())
  {
    const auto * node=entry.first;
    if(node->processor_id()!=processor_id())continue;
    const auto dof=node->dof_number(_sys.number(),_var->number(),0);
    ADReal pressure=(*_sys.currentSolution())(dof);Moose::derivInsert(pressure.derivatives(),dof,1.0);
    const ADReal residual=std::min(pressure,_c*entry.second.first/entry.second.second);
    addResidualsAndJacobian(_assembly,std::array<ADReal,1>{residual},std::array<dof_id_type,1>{dof},_var->scalingFactor());
  }
}
