#include "CDPVelocityDynamicPhysics.h"
#include "Factory.h"
#include "FEProblem.h"

registerMooseAction("DamSafetyApp", CDPVelocityDynamicPhysics, "meta_action");
registerMooseAction("DamSafetyApp", CDPVelocityDynamicPhysics, "setup_mesh_complete");
registerMooseAction("DamSafetyApp", CDPVelocityDynamicPhysics, "validate_coordinate_systems");
registerMooseAction("DamSafetyApp", CDPVelocityDynamicPhysics, "add_variable");
registerMooseAction("DamSafetyApp", CDPVelocityDynamicPhysics, "add_aux_variable");
registerMooseAction("DamSafetyApp", CDPVelocityDynamicPhysics, "add_kernel");
registerMooseAction("DamSafetyApp", CDPVelocityDynamicPhysics, "add_aux_kernel");
registerMooseAction("DamSafetyApp", CDPVelocityDynamicPhysics, "add_material");
registerMooseAction("DamSafetyApp", CDPVelocityDynamicPhysics, "add_master_action_material");

InputParameters
CDPVelocityDynamicPhysics::validParams()
{
  auto params = DynamicSolidMechanicsPhysics::validParams();
  params.addClassDescription("Default-off 3D diagnostic total-strain-rate damping; off delegates to upstream dynamics.");
  params.addParam<bool>("enable_cdp_velocity_damping", false, "Enable the diagnostic candidate");
  params.addParam<bool>("diagnostic_3d_model", false, "Explicit acknowledgement: 3D SMALL constant isotropic C0 and no eigenstrain");
  params.addParam<bool>("require_cdp_material", true, "False only for elastic diagnostic controls");
  params.addParam<Real>("expected_youngs_modulus", "Constant undamaged Young modulus, identical to the CDP integrator");
  params.addParam<Real>("expected_poissons_ratio", "Constant undamaged Poisson ratio, identical to the CDP integrator");
  params.addParam<Real>("expected_stiffness_damping", "Constant nonnegative stiffness damping coefficient");
  return params;
}

CDPVelocityDynamicPhysics::CDPVelocityDynamicPhysics(const InputParameters & parameters)
  : DynamicSolidMechanicsPhysics(parameters), _enabled(getParam<bool>("enable_cdp_velocity_damping"))
{
  if (_enabled && (!getParam<bool>("diagnostic_3d_model") || _ndisp != 3 ||
      getParam<MooseEnum>("strain") != "SMALL" ||
      !isParamValid("incremental") || !getParam<bool>("incremental") ||
      getParam<bool>("volumetric_locking_correction") || getParam<bool>("static_initialization") ||
      isParamValid("out_of_plane_strain") || _use_ad || _lagrangian_kernels ||
      !getParam<std::vector<MaterialPropertyName>>("eigenstrain_names").empty() ||
      !isParamValid("expected_youngs_modulus") || !isParamValid("expected_poissons_ratio") ||
      !isParamValid("expected_stiffness_damping")))
    mooseError("CDP velocity candidate requires acknowledged incremental SMALL 3D with constant C0, no AD/Bbar/thickness/static initialization");
}

std::string
CDPVelocityDynamicPhysics::getKernelType()
{
  const auto original = DynamicSolidMechanicsPhysics::getKernelType();
  return _enabled ? "CDPVelocityStressDivergence" : original;
}

InputParameters
CDPVelocityDynamicPhysics::getKernelParameters(std::string type)
{
  auto params = DynamicSolidMechanicsPhysics::getKernelParameters(type);
  if (_enabled && type == "CDPVelocityStressDivergence")
  {
    params.set<Real>("beta") = getParam<Real>("newmark_beta");
    params.set<Real>("gamma") = getParam<Real>("newmark_gamma");
  }
  return params;
}

void
CDPVelocityDynamicPhysics::act()
{
  DynamicSolidMechanicsPhysics::act();
  if (!_enabled || _current_task != "add_material")
    return;
  auto params = _factory.getValidParams("CDPVelocityRateMaterial");
  params.applyParameters(parameters(), {"displacements", "velocities", "accelerations"});
  params.set<std::vector<VariableName>>("displacements") = _displacements;
  params.set<std::vector<VariableName>>("velocities") = std::vector<VariableName>(_velocities.begin(), _velocities.begin() + _ndisp);
  params.set<std::vector<VariableName>>("accelerations") = std::vector<VariableName>(_accelerations.begin(), _accelerations.begin() + _ndisp);
  params.set<Real>("beta") = getParam<Real>("newmark_beta");
  params.set<Real>("gamma") = getParam<Real>("newmark_gamma");
  params.set<Real>("alpha") = getParam<Real>("hht_alpha");
  params.set<MaterialPropertyName>("zeta") = getParam<MaterialPropertyName>("stiffness_damping_coefficient");
  _problem->addMaterial("CDPVelocityRateMaterial", name() + "_cdp_velocity_material", params);
}
