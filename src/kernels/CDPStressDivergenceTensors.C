#include "CDPStressDivergenceTensors.h"
#include "CDPBbarContraction.h"
registerMooseObject("DamSafetyApp", CDPStressDivergenceTensors);
InputParameters CDPStressDivergenceTensors::validParams()
{
  auto p=StressDivergenceTensors::validParams();
  p.addClassDescription("Full B-bar Jacobian contraction for 3D small-strain CDP tangents.");
  return p;
}
CDPStressDivergenceTensors::CDPStressDivergenceTensors(const InputParameters & p)
  : StressDivergenceTensors(p)
{
  if (_ndisp != 3 || _use_finite_deform_jacobian || _use_displaced_mesh)
    mooseError("CDPStressDivergenceTensors requires 3D small strain on the undisplaced mesh");
}
Real CDPStressDivergenceTensors::bbarJacobian(unsigned int b) const
{
  const auto & g=_grad_test[_i][_qp];
  const auto & h=_grad_phi[_j][_qp];
  return CDPBbar::contract(_Jacobian_mult[_qp],_component,b,
                          {g(0),g(1),g(2)},{h(0),h(1),h(2)},
                          _avg_grad_test[_i][_component],_avg_grad_phi[_j][b]);
}
Real CDPStressDivergenceTensors::computeQpJacobian()
{
  return _volumetric_locking_correction ? bbarJacobian(_component)
                                       : StressDivergenceTensors::computeQpJacobian();
}
Real CDPStressDivergenceTensors::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_volumetric_locking_correction)
    for (unsigned int b=0;b<_ndisp;++b)
      if (jvar==_disp_var[b])
        return bbarJacobian(b);
  return StressDivergenceTensors::computeQpOffDiagJacobian(jvar);
}
