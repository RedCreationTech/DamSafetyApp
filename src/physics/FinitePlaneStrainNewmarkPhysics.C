#include "FinitePlaneStrainNewmarkPhysics.h"
#include "Factory.h"
#include "FEProblem.h"

registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "meta_action");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "setup_mesh_complete");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "validate_coordinate_systems");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_variable");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_aux_variable");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_kernel");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_aux_kernel");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_material");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_master_action_material");

InputParameters
FinitePlaneStrainNewmarkPhysics::validParams()
{
  auto params = DynamicSolidMechanicsPhysics::validParams();
  params.addClassDescription("Default-off finite plane-strain Newmark C0 damping; off delegates "
                             "exactly to upstream dynamics.");
  params.addParam<bool>("enable_finite_newmark_damping",
                        false,
                        "Replace only stiffness damping with the finite velocity-rate candidate");
  params.addParam<Real>("expected_youngs_modulus", "Constant C0 Young modulus");
  params.addParam<Real>("expected_poissons_ratio", "Constant C0 Poisson ratio");
  params.addParam<Real>("expected_stiffness_damping",
                        "Constant nonnegative original stiffness damping");
  return params;
}

FinitePlaneStrainNewmarkPhysics::FinitePlaneStrainNewmarkPhysics(const InputParameters & parameters)
  : DynamicSolidMechanicsPhysics(parameters),
    _enabled(getParam<bool>("enable_finite_newmark_damping"))
{
  if (_enabled &&
      (_ndisp != 2 || !_use_displaced_mesh || getParam<MooseEnum>("strain") != "FINITE" ||
       getParam<MooseEnum>("planar_formulation") != "PLANE_STRAIN" ||
       isParamValid("out_of_plane_strain") || _use_ad || _lagrangian_kernels ||
       getParam<bool>("volumetric_locking_correction") || getParam<bool>("static_initialization") ||
       isParamValid("base_name") ||
       !getParam<std::vector<MaterialPropertyName>>("eigenstrain_names").empty() ||
       !isParamValid("expected_youngs_modulus") || !isParamValid("expected_poissons_ratio") ||
       !isParamValid("expected_stiffness_damping")))
    mooseError("Finite Newmark candidate requires FINITE PLANE_STRAIN, constant C0, standard "
               "non-AD stress kernels, no Bbar/static initialization/base name/eigenstrain");
}

InputParameters
FinitePlaneStrainNewmarkPhysics::getKernelParameters(std::string type)
{
  auto params = DynamicSolidMechanicsPhysics::getKernelParameters(type);
  if (_enabled && type == "DynamicStressDivergenceTensors")
    params.set<MaterialPropertyName>("zeta") = "0.0";
  return params;
}

void
FinitePlaneStrainNewmarkPhysics::act()
{
  DynamicSolidMechanicsPhysics::act();
  if (!_enabled || _current_task != "add_kernel")
    return;
  const std::array<std::string, 2> axes{{"x", "y"}};
  for (unsigned int c = 0; c < 2; ++c)
  {
    auto params = _factory.getValidParams("ADFinitePlaneStrainNewmarkDamping");
    params.applySpecificParameters(parameters(),
                                   {"block",
                                    "expected_youngs_modulus",
                                    "expected_poissons_ratio",
                                    "expected_stiffness_damping"});
    params.set<NonlinearVariableName>("variable") = _displacements[c];
    params.set<unsigned int>("component") = c;
    params.set<std::vector<VariableName>>("displacements") = _displacements;
    params.set<std::vector<VariableName>>("velocities") =
        std::vector<VariableName>(_velocities.begin(), _velocities.begin() + 2);
    params.set<std::vector<VariableName>>("accelerations") =
        std::vector<VariableName>(_accelerations.begin(), _accelerations.begin() + 2);
    params.set<Real>("beta") = getParam<Real>("newmark_beta");
    params.set<Real>("gamma") = getParam<Real>("newmark_gamma");
    params.set<Real>("alpha") = getParam<Real>("hht_alpha");
    params.set<MaterialPropertyName>("zeta") =
        getParam<MaterialPropertyName>("stiffness_damping_coefficient");
    _problem->addKernel("ADFinitePlaneStrainNewmarkDamping",
                        "ds04g_" + name() + "_velocity_damping_" + axes[c],
                        params);
  }
}
