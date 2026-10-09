#include "ElasticVelocityDynamicPhysics.h"
#include "Factory.h"
#include "FEProblem.h"

registerMooseAction("DamSafetyApp", ElasticVelocityDynamicPhysics, "meta_action");
registerMooseAction("DamSafetyApp", ElasticVelocityDynamicPhysics, "setup_mesh_complete");
registerMooseAction("DamSafetyApp", ElasticVelocityDynamicPhysics, "validate_coordinate_systems");
registerMooseAction("DamSafetyApp", ElasticVelocityDynamicPhysics, "add_variable");
registerMooseAction("DamSafetyApp", ElasticVelocityDynamicPhysics, "add_aux_variable");
registerMooseAction("DamSafetyApp", ElasticVelocityDynamicPhysics, "add_kernel");
registerMooseAction("DamSafetyApp", ElasticVelocityDynamicPhysics, "add_aux_kernel");
registerMooseAction("DamSafetyApp", ElasticVelocityDynamicPhysics, "add_material");
registerMooseAction("DamSafetyApp", ElasticVelocityDynamicPhysics, "add_master_action_material");

InputParameters
ElasticVelocityDynamicPhysics::validParams()
{
  auto params = DynamicSolidMechanicsPhysics::validParams();
  params.addClassDescription("Default-off strictly elastic Newmark velocity damping candidate. Off uses the original upstream kernels and adds no rate equation.");
  params.addParam<bool>("enable_velocity_rate_damping", false, "Enable elastic-only velocity damping candidate");
  params.addParam<bool>("strict_elastic_model", false, "Explicitly acknowledge the audited constant linear elastic material and no eigenstrain/damage");
  params.addParam<VariableName>("thickness_rate", "Existing nonlinear thickness rate variable in the original thickness FE/block space");
  params.addParam<MooseEnum>("thickness_rate_formulation", MooseEnum("same_space pointwise", "same_space"),
                           "Same-space algebraic variable or independently audited observable elimination");
  params.addParam<bool>("pointwise_rate_equivalence_audited", false,
                       "Explicitly acknowledge full-row-rank observable plane-stress rate closure for this input");
  params.addParam<Real>("expected_youngs_modulus", "Audited constant Young modulus");
  params.addParam<Real>("expected_poissons_ratio", "Audited constant Poisson ratio");
  return params;
}

ElasticVelocityDynamicPhysics::ElasticVelocityDynamicPhysics(const InputParameters & parameters)
  : DynamicSolidMechanicsPhysics(parameters), _enabled(getParam<bool>("enable_velocity_rate_damping"))
{
  if (_enabled)
  {
    if (!getParam<bool>("strict_elastic_model") || getParam<MooseEnum>("strain") != "SMALL" ||
        getParam<MooseEnum>("planar_formulation") != "WEAK_PLANE_STRESS" ||
        !isParamValid("incremental") || !getParam<bool>("incremental") ||
        getParam<bool>("volumetric_locking_correction") || getParam<bool>("static_initialization") ||
        getParam<MooseEnum>("out_of_plane_direction") != "z" ||
        !getParam<std::vector<MaterialPropertyName>>("eigenstrain_names").empty() ||
        _use_ad || _lagrangian_kernels || _ndisp != 2 ||
        !isParamValid("expected_youngs_modulus") ||
        !isParamValid("expected_poissons_ratio"))
      mooseError("Elastic velocity candidate requires acknowledged incremental SMALL WEAK_PLANE_STRESS elasticity, explicit rate variable and elastic constants, no AD/Bbar/static initialization");
    const bool pointwise = getParam<MooseEnum>("thickness_rate_formulation") == "pointwise";
    if ((pointwise && (!getParam<bool>("pointwise_rate_equivalence_audited") ||
                      isParamValid("thickness_rate"))) ||
        (!pointwise && !isParamValid("thickness_rate")))
      mooseError("Pointwise rate requires explicit audited equivalence and no thickness_rate variable; same_space requires thickness_rate");
  }
}

std::string
ElasticVelocityDynamicPhysics::getKernelType()
{
  const auto original = DynamicSolidMechanicsPhysics::getKernelType();
  return _enabled ? "ElasticVelocityStressDivergence" : original;
}

InputParameters
ElasticVelocityDynamicPhysics::getKernelParameters(std::string type)
{
  auto params = DynamicSolidMechanicsPhysics::getKernelParameters(type);
  if (_enabled && (type == "ElasticVelocityStressDivergence" || type == "ElasticPlaneStressRate"))
  {
    params.set<Real>("beta") = getParam<Real>("newmark_beta");
    params.set<Real>("gamma") = getParam<Real>("newmark_gamma");
    if (type == "ElasticVelocityStressDivergence" &&
        getParam<MooseEnum>("thickness_rate_formulation") == "same_space")
      params.set<std::vector<VariableName>>("out_of_plane_strain_rate") = {getParam<VariableName>("thickness_rate")};
  }
  return params;
}

void
ElasticVelocityDynamicPhysics::act()
{
  DynamicSolidMechanicsPhysics::act();
  if (!_enabled)
    return;
  if (_current_task == "add_kernel" &&
      getParam<MooseEnum>("thickness_rate_formulation") == "same_space")
  {
    auto params = getKernelParameters("ElasticPlaneStressRate");
    params.set<NonlinearVariableName>("variable") = getParam<VariableName>("thickness_rate");
    _problem->addKernel("ElasticPlaneStressRate", name() + "_elastic_rate_constraint", params);
  }
  if (_current_task == "add_material")
  {
    auto params = _factory.getValidParams("ElasticVelocityRateMaterial");
    params.applyParameters(parameters(), {"displacements", "velocities", "accelerations", "out_of_plane_strain"});
    params.set<std::vector<VariableName>>("displacements") = _displacements;
    params.set<std::vector<VariableName>>("velocities") = std::vector<VariableName>(_velocities.begin(), _velocities.begin() + _ndisp);
    params.set<std::vector<VariableName>>("accelerations") = std::vector<VariableName>(_accelerations.begin(), _accelerations.begin() + _ndisp);
    params.set<std::vector<VariableName>>("out_of_plane_strain") = {getParam<VariableName>("out_of_plane_strain")};
    if (getParam<MooseEnum>("thickness_rate_formulation") == "same_space")
      params.set<std::vector<VariableName>>("out_of_plane_strain_rate") = {getParam<VariableName>("thickness_rate")};
    params.set<bool>("pointwise_plane_stress_rate") =
        getParam<MooseEnum>("thickness_rate_formulation") == "pointwise";
    params.set<Real>("beta") = getParam<Real>("newmark_beta");
    params.set<Real>("gamma") = getParam<Real>("newmark_gamma");
    params.set<Real>("alpha") = getParam<Real>("hht_alpha");
    params.set<MaterialPropertyName>("zeta") = getParam<MaterialPropertyName>("stiffness_damping_coefficient");
    _problem->addMaterial("ElasticVelocityRateMaterial", name() + "_elastic_rate_material", params);
  }
}
