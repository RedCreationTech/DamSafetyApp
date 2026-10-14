#include "CDPPlaneStressShearJacobian.h"
#include "PlaneStressCoupling.h"
#include "RankFourTensor.h"
#include "MooseMesh.h"

registerMooseObject("DamSafetyApp", CDPPlaneStressShearJacobian);

InputParameters
CDPPlaneStressShearJacobian::validParams()
{
  auto params = WeakPlaneStress::validParams();
  params.addClassDescription("Default-off Jacobian-only XY shear completion; zero residual and diagonal, no HHT factor on the static thickness row.");
  params.addParam<bool>("complete_plane_stress_coupling", false, "Enable only the missing displacement columns");
  return params;
}

CDPPlaneStressShearJacobian::CDPPlaneStressShearJacobian(const InputParameters & parameters)
  : WeakPlaneStress(parameters),
    _complete_coupling(getParam<bool>("complete_plane_stress_coupling"))
{
  if (_complete_coupling && (_mesh.dimension() != 2 || _ndisp != 2 || _direction != 2 ||
      _volumetric_locking_correction || _temp_coupled || getParam<bool>("use_displaced_mesh") ||
      !getParam<std::vector<MaterialPropertyName>>("eigenstrain_names").empty() ||
      _var.feType().family != LAGRANGE || _var.feType().order != FIRST))
    mooseError("Complete plane-stress coupling requires XY small-strain FIRST LAGRANGE thickness without Bbar, temperature or eigenstrain");
}

Real
CDPPlaneStressShearJacobian::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (!_complete_coupling)
    return 0.;
  for (unsigned int component = 0; component < 2; ++component)
    if (jvar == _disp_var[component])
      return PlaneStressCoupling::thicknessShearRow(
          _Jacobian_mult[_qp], component, _test[_i][_qp], _grad_phi[_j][_qp]);
  return 0.;
}
