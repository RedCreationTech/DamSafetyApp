#include "AqueductHHTInertia.h"
#include "AqueductHHTScheme.h"
registerMooseObject("DamSafetyApp",AqueductHHTInertia);
InputParameters AqueductHHTInertia::validParams()
{
  auto p=InertialForce::validParams();
  p.addRequiredRangeCheckedParam<Real>("hht_alpha","hht_alpha>=-0.5 & hht_alpha<=0","Same whole-system HHT alpha");
  p.addParam<Real>("hht_start_time",0.,"No inertia during source static stage");
  p.suppressParameter<Real>("alpha");p.suppressParameter<MaterialPropertyName>("eta");
  return p;
}
AqueductHHTInertia::AqueductHHTInertia(const InputParameters & p)
 :InertialForce(p),_hht_mass_scale(AqueductHHTScheme::massScale(getParam<Real>("hht_alpha"))),
 _hht_start_time(getParam<Real>("hht_start_time")){}
Real AqueductHHTInertia::computeQpResidual()
{return _t<=_hht_start_time?0.:_hht_mass_scale*InertialForce::computeQpResidual();}
Real AqueductHHTInertia::computeQpJacobian()
{return _t<=_hht_start_time?0.:_hht_mass_scale*InertialForce::computeQpJacobian();}
