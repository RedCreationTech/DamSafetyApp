#include "AqueductPhaseDT.h"
#include <algorithm>
registerMooseObject("DamSafetyApp",AqueductPhaseDT);
InputParameters AqueductPhaseDT::validParams()
{
  auto p=FunctionDT::validParams();
  p.addRequiredRangeCheckedParam<Real>("static_end_time","static_end_time>0","Accepted static endpoint");
  for(const auto * n:{"static_dtmin","static_dtmax","dynamic_dtmin","dynamic_dtmax"})
    p.addRequiredRangeCheckedParam<Real>(n,std::string(n)+">0","Source phase timestep bound");
  return p;
}
AqueductPhaseDT::AqueductPhaseDT(const InputParameters & p)
 :FunctionDT(p),_static_end(getParam<Real>("static_end_time")),
 _static_min(getParam<Real>("static_dtmin")),_static_max(getParam<Real>("static_dtmax")),
 _dynamic_min(getParam<Real>("dynamic_dtmin")),_dynamic_max(getParam<Real>("dynamic_dtmax"))
{
  if(_static_min>_static_max || _dynamic_min>_dynamic_max ||
     _dt_min>std::min(_static_min,_dynamic_min) || _dt_max<std::max(_static_max,_dynamic_max))
    paramError("static_dtmin","Phase bounds must be ordered and contained in Executioner dtmin/dtmax");
}
bool AqueductPhaseDT::constrainStep(Real & dt)
{
  const bool dynamic=_time>=_static_end;
  dt=std::max(dynamic?_dynamic_min:_static_min,std::min(dt,dynamic?_dynamic_max:_static_max));
  // The parent may reduce the final increment to a physical endpoint.
  return FunctionDT::constrainStep(dt);
}
Real AqueductPhaseDT::computeFailedDT()
{
  const Real minimum=_time_old>=_static_end?_dynamic_min:_static_min;
  if(_dt<=minimum)mooseError("Aqueduct phase solve failed at source dtmin; no tolerance relaxation");
  return std::max(minimum,_cutback_factor_at_failure*_dt);
}
