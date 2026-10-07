#include "GroundMotionInertia.h"
#include "Function.h"
registerMooseObject("DamSafetyApp", GroundMotionInertia);
InputParameters GroundMotionInertia::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription("Uniform prescribed translation in relative-coordinate inertia, including HHT-weighted mass damping.");
  params.addRequiredParam<FunctionName>("ground_acceleration", "Analytic ground acceleration");
  params.addRequiredParam<FunctionName>("ground_velocity", "Analytic ground velocity");
  params.addParam<MaterialPropertyName>("density", "density", "Density");
  params.addParam<MaterialPropertyName>("eta", 0.0, "Mass Rayleigh damping coefficient");
  params.addParam<Real>("alpha", 0.0, "HHT alpha");
  params.set<bool>("use_displaced_mesh") = false;
  return params;
}
GroundMotionInertia::GroundMotionInertia(const InputParameters & parameters)
  : Kernel(parameters), _acceleration(getFunction("ground_acceleration")),
    _velocity(getFunction("ground_velocity")), _density(getMaterialProperty<Real>("density")),
    _eta(getMaterialProperty<Real>("eta")), _alpha(getParam<Real>("alpha")) {}
Real GroundMotionInertia::computeQpResidual()
{
  if (_dt <= 0) return 0;
  return _test[_i][_qp] * _density[_qp] *
         (_acceleration.value(_t, _q_point[_qp]) + _eta[_qp] *
          ((1+_alpha)*_velocity.value(_t, _q_point[_qp]) -
           _alpha*_velocity.value(_t-_dt, _q_point[_qp])));
}
