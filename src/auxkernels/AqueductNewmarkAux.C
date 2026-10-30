#include "AqueductNewmarkAux.h"

registerMooseObject("DamSafetyApp",AqueductNewmarkAcceleration);
registerMooseObject("DamSafetyApp",AqueductNewmarkVelocity);

InputParameters AqueductNewmarkAcceleration::validParams()
{
  auto p=NewmarkAccelAux::validParams();
  p.addRequiredRangeCheckedParam<Real>("static_end_time","static_end_time>=0",
                                      "Same accepted static endpoint as AqueductHHT");
  return p;
}
AqueductNewmarkAcceleration::AqueductNewmarkAcceleration(const InputParameters & p)
  :NewmarkAccelAux(p),_static_end(getParam<Real>("static_end_time")) {}
Real AqueductNewmarkAcceleration::computeValue()
{
  return _t<=_static_end?0.:NewmarkAccelAux::computeValue();
}

InputParameters AqueductNewmarkVelocity::validParams()
{
  auto p=NewmarkVelAux::validParams();
  p.addRequiredRangeCheckedParam<Real>("static_end_time","static_end_time>=0",
                                      "Same accepted static endpoint as AqueductHHT");
  return p;
}
AqueductNewmarkVelocity::AqueductNewmarkVelocity(const InputParameters & p)
  :NewmarkVelAux(p),_static_end(getParam<Real>("static_end_time")) {}
Real AqueductNewmarkVelocity::computeValue()
{
  return _t<=_static_end?0.:NewmarkVelAux::computeValue();
}
