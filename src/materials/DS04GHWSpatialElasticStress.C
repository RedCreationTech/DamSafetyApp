#include "DS04GHWSpatialElasticStress.h"

registerMooseObject("DamSafetyApp", DS04GHWSpatialElasticStress);

InputParameters
DS04GHWSpatialElasticStress::validParams()
{
  InputParameters params = ComputeFiniteStrainElasticStress::validParams();
  params.addParam<bool>("acknowledge_research_limits", false,
                        "Explicitly accept elastic-only diagnostic behavior and approximate tangent");
  params.addClassDescription("Research only: isotropic elastic spatial Hughes-Winget update; "
                             "no CDP, no exact finite-strain tangent.");
  return params;
}

DS04GHWSpatialElasticStress::DS04GHWSpatialElasticStress(const InputParameters & parameters)
  : ComputeFiniteStrainElasticStress(parameters),
    _hw_midpoint_gradient(getMaterialPropertyByName<RankTwoTensor>(_base_name + "def_grad_mid")),
    _hw_displacement_increment_gradient(getMaterialPropertyByName<RankTwoTensor>(_base_name + "f_bar"))
{
  if (!getParam<bool>("acknowledge_research_limits"))
    paramError("acknowledge_research_limits", "This elastic-only diagnostic is not a production material");
}

void
DS04GHWSpatialElasticStress::initialSetup()
{
  ComputeFiniteStrainElasticStress::initialSetup();
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("DS04GHWSpatialElasticStress is restricted to guaranteed isotropic elasticity");
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::CONSTANT_IN_TIME))
    mooseError("DS04GHWSpatialElasticStress requires time-constant elasticity");
  if (hasMaterialPropertyByName<RankTwoTensor>(_base_name + "extra_stress"))
    mooseError("DS04GHWSpatialElasticStress does not support extra_stress");
}

void
DS04GHWSpatialElasticStress::computeQpStress()
{
  // The two named properties are dependencies. Their presence is not a proof of
  // producer identity; users must independently lock Hughes-Winget kinematics.
  (void)_hw_midpoint_gradient;
  (void)_hw_displacement_increment_gradient;
  const auto & rotation = _rotation_increment[_qp];
  _stress[_qp] = rotation * _stress_old[_qp] * rotation.transpose() +
                 _elasticity_tensor[_qp] * _strain_increment[_qp];
  _elastic_strain[_qp] = rotation * _elastic_strain_old[_qp] * rotation.transpose() +
                         _strain_increment[_qp];
  _rotation_total[_qp] = rotation * _rotation_total_old[_qp];
  // Elastic approximation only; not the exact finite-deformation Jacobian.
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}
