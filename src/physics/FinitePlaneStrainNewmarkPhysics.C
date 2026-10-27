#include "ActionWarehouse.h"
#include "FEProblem.h"
#include "Factory.h"
#include "FinitePlaneStrainNewmarkPhysics.h"
#include "MooseMesh.h"
#include "MooseObjectAction.h"
#include "SetupQuadratureAction.h"
#include "libmesh/enum_quadrature_type.h"
#include <cmath>

registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "meta_action");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "setup_mesh_complete");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "validate_coordinate_systems");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_variable");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_aux_variable");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_kernel");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_aux_kernel");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_material");
registerMooseAction("DamSafetyApp", FinitePlaneStrainNewmarkPhysics, "add_master_action_material");

InputParameters FinitePlaneStrainNewmarkPhysics::validParams()
{
  auto params = DynamicSolidMechanicsPhysics::validParams();
  params.addClassDescription("Default-off finite plane-strain Newmark C0 damping; off delegates "
                             "exactly to upstream dynamics.");
  params.addParam<bool>("enable_finite_newmark_damping", false,
                        "Replace only stiffness damping with the finite velocity-rate candidate");
  params.addParam<Real>("expected_youngs_modulus", "Constant C0 Young modulus");
  params.addParam<Real>("expected_poissons_ratio", "Constant C0 Poisson ratio");
  params.addParam<Real>("expected_stiffness_damping",
                        "Constant nonnegative original stiffness damping");
  params.addParam<MooseEnum>(
      "spatial_discretization", MooseEnum("legacy_q4 sp1_center_hg", "legacy_q4"),
      "Opt-in DS04G center entity plus reference HG; never a CDP formula change");
  params.addParam<Real>("sp1_density", "Required frozen constant density for SP1");
  params.addParam<Real>("sp1_static_end", 1.000001,
                        "Must match the TimePeriod static phase endpoint");
  params.addParam<Real>("sp1_gravity", -9.8,
                        "Uniform Y gravity; must match original Gravity kernel");
  params.addParam<bool>("sp1_consistent_tangent", false,
                        "Opt-in rotation-false Taylor/geometry Jacobian; original residual retained");
  return params;
}

FinitePlaneStrainNewmarkPhysics::FinitePlaneStrainNewmarkPhysics(const InputParameters &parameters)
    : DynamicSolidMechanicsPhysics(parameters),
      _enabled(getParam<bool>("enable_finite_newmark_damping")),
      _sp1(getParam<MooseEnum>("spatial_discretization") == "sp1_center_hg"),
      _sp1_consistent(getParam<bool>("sp1_consistent_tangent"))
{
  if (_sp1_consistent && !_sp1)
    paramError("sp1_consistent_tangent", "Requires spatial_discretization=sp1_center_hg");
  if (_sp1_consistent && isParamValid("decomposition_method") &&
      getParam<MooseEnum>("decomposition_method") != "TaylorExpansion")
    paramError("sp1_consistent_tangent", "Requires the frozen TaylorExpansion strain path");
  if (_sp1 &&
      (_enabled || _ndisp != 2 || !_use_displaced_mesh ||
       getParam<MooseEnum>("strain") != "FINITE" ||
       getParam<MooseEnum>("planar_formulation") != "PLANE_STRAIN" ||
       isParamValid("out_of_plane_strain") || _use_ad || _lagrangian_kernels ||
       getParam<bool>("volumetric_locking_correction") || getParam<bool>("static_initialization") ||
       isParamValid("base_name") ||
       !getParam<std::vector<MaterialPropertyName>>("eigenstrain_names").empty() ||
       !isParamValid("expected_youngs_modulus") || !isParamValid("expected_poissons_ratio") ||
       !isParamValid("expected_stiffness_damping") || !isParamValid("sp1_density")))
    mooseError("SP1 requires FINITE PLANE_STRAIN, constant E/nu/density, Newmark damping OFF, no "
               "Bbar/extra strain/eigenstrain/base_name");
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

InputParameters FinitePlaneStrainNewmarkPhysics::getKernelParameters(std::string type)
{
  auto params = DynamicSolidMechanicsPhysics::getKernelParameters(type);
  if (_enabled && type == "DynamicStressDivergenceTensors")
    params.set<MaterialPropertyName>("zeta") = "0.0";
  if (_sp1_consistent && type == "DS04GSP1StressDivergence")
    params.set<bool>("sp1_dynamic") = true;
  return params;
}

std::string FinitePlaneStrainNewmarkPhysics::getKernelType()
{
  return _sp1_consistent ? "DS04GSP1StressDivergence"
                         : DynamicSolidMechanicsPhysics::getKernelType();
}

void FinitePlaneStrainNewmarkPhysics::act()
{
  // Set the volume rule after the standard setup_quadrature task but before
  // material state initialization. The legacy branch never calls createQRules.
  if (_sp1 && _current_task == "add_material")
  {
    if (_mesh->dimension() != 2 || _mesh->meshSubdomains().size() != 1)
      mooseError("SP1 admission is a single-subdomain 2D QUAD4 mesh");
    for (const auto *q : _awh.getActions<SetupQuadratureAction>())
      for (const auto &n : {"type", "order", "element_order", "side_order", "custom_blocks",
                            "custom_orders", "custom_types"})
        if (q->parameters().isParamSetByUser(n))
          mooseError(
              "SP1 owns the center volume quadrature; do not specify a separate Quadrature block");
    _problem->createQRules(libMesh::QGAUSS, libMesh::INVALID_ORDER, libMesh::FIRST,
                           libMesh::INVALID_ORDER);
  }
  DynamicSolidMechanicsPhysics::act();
  if (_sp1 && _current_task == "add_kernel")
  {
    auto constant = [&](const std::string &name)
    {
      const auto text = getParam<MaterialPropertyName>(name);
      std::size_t used = 0;
      double value = 0;
      try
      {
        value = std::stod(text, &used);
      }
      catch (...)
      {
        mooseError("SP1 requires numeric constant ", name);
      }
      if (used != text.size() || !std::isfinite(value))
        mooseError("SP1 requires finite numeric ", name);
      return value;
    };
    unsigned gravity_count = 0, phase_count = 0;
    for (const auto *action : _awh.getActions<MooseObjectAction>())
    {
      const auto &q = action->getObjectParams();
      if (action->getMooseObjectType() == "Gravity")
      {
        ++gravity_count;
        if (q.get<NonlinearVariableName>("variable") != _displacements[1] ||
            std::abs(q.get<Real>("value") - getParam<Real>("sp1_gravity")) > 1e-14 ||
            std::abs(q.get<Real>("alpha") - getParam<Real>("hht_alpha")) > 1e-14 ||
            q.get<FunctionName>("function") != "1" || !q.get<bool>("use_displaced_mesh"))
          mooseError("SP1 requires the original uniform current-geometry Y Gravity and matching "
                     "alpha/value");
      }
      if (action->getMooseObjectType() == "TimePeriod")
      {
        ++phase_count;
        if (!q.isParamValid("start_time") || !q.isParamValid("end_time") ||
            q.get<std::vector<Real>>("start_time") != std::vector<Real>{0.} ||
            q.get<std::vector<Real>>("end_time") !=
                std::vector<Real>{getParam<Real>("sp1_static_end")})
          mooseError("SP1 static_end must match the original single TimePeriod phase");
      }
      if (action->getMooseObjectType() == "ComputeMultipleInelasticStress" &&
          q.get<bool>("perform_finite_strain_rotations"))
        mooseError(
            "SP1 retains original finite-strain rotation policy; do not change it simultaneously");
    }
    if (gravity_count != 1 || phase_count != 1)
      mooseError("SP1 requires exactly one original uniform Gravity and one TimePeriod phase");
    const double eta = constant("mass_damping_coefficient");
    const double zeta = constant("stiffness_damping_coefficient");
    if (std::abs(zeta - getParam<Real>("expected_stiffness_damping")) > 1e-14)
      mooseError("SP1 declared stiffness damping mismatch");
    for (unsigned c = 0; c < 2; ++c)
      for (const auto &term : {"hourglass", "mass", "gravity"})
      {
        if (std::string(term) == "gravity" && c == 0)
          continue;
        auto p = _factory.getValidParams("ADDS04GSP1Correction");
        p.applySpecificParameters(parameters(), {"block"});
        p.set<NonlinearVariableName>("variable") = _displacements[c];
        p.set<std::vector<VariableName>>("displacements") = _displacements;
        p.set<std::vector<VariableName>>("velocity") = {_velocities[c]};
        p.set<std::vector<VariableName>>("acceleration") = {_accelerations[c]};
        p.set<MooseEnum>("term") = term;
        p.set<unsigned>("component") = c;
        p.set<Real>("youngs_modulus") = getParam<Real>("expected_youngs_modulus");
        p.set<Real>("poissons_ratio") = getParam<Real>("expected_poissons_ratio");
        p.set<Real>("expected_density") = getParam<Real>("sp1_density");
        p.set<Real>("eta") = eta;
        p.set<Real>("zeta") = zeta;
        p.set<Real>("alpha") = getParam<Real>("hht_alpha");
        p.set<Real>("beta") = getParam<Real>("newmark_beta");
        p.set<Real>("gamma") = getParam<Real>("newmark_gamma");
        p.set<Real>("static_end") = getParam<Real>("sp1_static_end");
        p.set<Real>("gravity") = getParam<Real>("sp1_gravity");
        _problem->addKernel("ADDS04GSP1Correction",
                            "ds04g_sp1_" + std::string(term) + "_" + std::to_string(c), p);
      }
  }
  if (!_enabled || _current_task != "add_kernel")
    return;
  const std::array<std::string, 2> axes{{"x", "y"}};
  for (unsigned int c = 0; c < 2; ++c)
  {
    auto params = _factory.getValidParams("ADFinitePlaneStrainNewmarkDamping");
    params.applySpecificParameters(parameters(),
                                   {"block", "expected_youngs_modulus", "expected_poissons_ratio",
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
                        "ds04g_" + name() + "_velocity_damping_" + axes[c], params);
  }
}
