#include "DS04GSP1StressDivergence.h"
#include "DS04GFalseRotationTangent.h"
#include "ActionWarehouse.h"
#include "MooseApp.h"
#include "MooseMesh.h"
#include "MooseObjectAction.h"
#include "libmesh/elem.h"

registerMooseObject("DamSafetyApp", DS04GSP1StressDivergence);
InputParameters DS04GSP1StressDivergence::validParams()
{
  auto p = DynamicStressDivergenceTensors::validParams();
  p.addClassDescription("Opt-in DS04G center-QP Taylor/rotation-false compatible Jacobian; "
                        "original residual unchanged");
  p.addParam<bool>("sp1_dynamic", false, "Use the original dynamic stress residual and weights");
  return p;
}
DS04GSP1StressDivergence::DS04GSP1StressDivergence(const InputParameters & p)
    : DynamicStressDivergenceTensors(p), _sp1_dynamic(getParam<bool>("sp1_dynamic")),
      _sp1_F(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
      _sp1_Fold(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient"))
{
  if (_ndisp != 2 || !_use_displaced_mesh || _use_finite_deform_jacobian ||
      _volumetric_locking_correction || _out_of_plane_strain_coupled ||
      !_base_name.empty() || _static_initialization || _component > 1)
    mooseError("SP1 consistent tangent requires XY displaced FINITE, no Bbar/base_name/standard "
               "finite Jacobian/static initialization/out-of-plane coupling");
}
void DS04GSP1StressDivergence::initialSetup()
{
  DynamicStressDivergenceTensors::initialSetup();
  unsigned stress_models = 0;
  for (const auto * a : _app.actionWarehouse().getActions<MooseObjectAction>())
    if (a->getMooseObjectType() == "ComputeMultipleInelasticStress")
    {
      ++stress_models;
      if (a->getObjectParams().get<bool>("perform_finite_strain_rotations"))
        mooseError("SP1 consistent tangent requires perform_finite_strain_rotations=false");
    }
  if (stress_models != 1 || _mesh.dimension() != 2 || getBlockCoordSystem() != Moose::COORD_XYZ)
    mooseError("SP1 consistent tangent requires one rotation-false stress material in Cartesian XY");
}
Real DS04GSP1StressDivergence::computeQpResidual()
{
  return _sp1_dynamic ? DynamicStressDivergenceTensors::computeQpResidual()
                      : StressDivergenceTensors::computeQpResidual();
}
Real DS04GSP1StressDivergence::tangent(unsigned direction)
{
  if (_current_elem->type() != QUAD4 || _qrule->n_points() != 1)
    mooseError("SP1 consistent tangent requires one center QP on QUAD4");
  Real A = 1;
  RankTwoTensor effective = _stress[_qp];
  if (_sp1_dynamic)
  {
    if (_dt <= 0)
      return 0;
    A = 1 + _alpha + (1 + _alpha) * _zeta[_qp] / _dt;
    const Real B = -_alpha - (1 + 2 * _alpha) * _zeta[_qp] / _dt;
    const Real C = _alpha * _zeta[_qp] / _dt;
    effective = A * _stress[_qp] + B * _stress_old[_qp] + C * _stress_older[_qp];
  }
  const auto p = DS04GFalseRotation::contract(_sp1_F[_qp], _sp1_Fold[_qp], _Jacobian_mult[_qp],
      effective, A, _component, direction, _grad_test[_i][_qp], _grad_phi[_j][_qp]);
  return p.material + p.geometry;
}
Real DS04GSP1StressDivergence::computeQpJacobian() { return tangent(_component); }
Real DS04GSP1StressDivergence::computeQpOffDiagJacobian(unsigned jvar)
{
  for (unsigned b = 0; b < _ndisp; ++b)
    if (jvar == _disp_var[b])
      return tangent(b);
  return 0;
}
