#include "AbaqusCDPStressUpdate.h"

#include "MooseException.h"

#include <array>
#include <exception>

registerMooseObject("DamSafetyApp", AbaqusCDPStressUpdate);

InputParameters
AbaqusCDPStressUpdate::validParams()
{
  InputParameters params = StressUpdateBase::validParams();
  params.addClassDescription(
      "Table-driven small-strain Abaqus-CDP-compatible stress update with transactional "
      "substeps and an algorithmic tangent.");
  params.addRequiredParam<FileName>("compression_hardening_file", "Compression hardening CSV");
  params.addRequiredParam<FileName>("compression_damage_file", "Compression damage CSV");
  params.addRequiredParam<FileName>("tension_stiffening_file", "Tension stiffening CSV");
  params.addRequiredParam<FileName>("tension_damage_file", "Tension damage CSV");
  params.addRequiredRangeCheckedParam<Real>(
      "youngs_modulus", "youngs_modulus > 0", "Young's modulus used by the CDP integrator");
  params.addRequiredRangeCheckedParam<Real>(
      "poissons_ratio", "poissons_ratio > -1 & poissons_ratio < 0.5", "Poisson's ratio");
  params.addParam<Real>("dilation_angle", 36.31, "Dilation angle in degrees");
  params.addParam<Real>("eccentricity", 0.1, "Hyperbolic flow-potential eccentricity");
  params.addParam<Real>("biaxial_to_uniaxial_compression_ratio", 1.16, "fb0/fc0");
  params.addParam<Real>("tensile_meridian_ratio", 0.667, "Kc tensile meridian ratio");
  params.addRangeCheckedParam<Real>(
      "tension_recovery", 1.0, "tension_recovery >= 0 & tension_recovery <= 1", "w_t");
  params.addRangeCheckedParam<Real>(
      "compression_recovery", 0.0, "compression_recovery >= 0 & compression_recovery <= 1", "w_c");
  params.addRangeCheckedParam<Real>(
      "viscosity", 0.0, "viscosity >= 0", "Duvaut-Lions relaxation time");
  params.addRangeCheckedParam<unsigned int>(
      "maximum_local_iterations", 40, "maximum_local_iterations > 0", "Local Newton limit");
  params.addRangeCheckedParam<Real>(
      "local_residual_tolerance", 1.0e-9, "local_residual_tolerance > 0", "Local tolerance");
  params.addRangeCheckedParam<Real>("local_finite_difference_step",
                                    1.0e-7,
                                    "local_finite_difference_step > 0",
                                    "Local Jacobian perturbation");
  params.addRangeCheckedParam<Real>("minimum_line_search",
                                    1.0e-6,
                                    "minimum_line_search > 0 & minimum_line_search <= 1",
                                    "Minimum local line-search factor");
  params.addRangeCheckedParam<unsigned int>(
      "maximum_substeps", 256, "maximum_substeps > 0", "Maximum binary material substeps");
  params.addRangeCheckedParam<Real>("maximum_strain_increment",
                                    0.0,
                                    "maximum_strain_increment >= 0",
                                    "Optional proactive material strain-increment limit");
  params.addRangeCheckedParam<Real>("reference_tangent_perturbation",
                                    1.0e-8,
                                    "reference_tangent_perturbation > 0",
                                    "Diagnostic reference-tangent perturbation");
  return params;
}

AbaqusCDPStressUpdate::AbaqusCDPStressUpdate(const InputParameters & parameters)
  : StressUpdateBase(parameters),
    _table(getParam<FileName>("compression_hardening_file"),
           getParam<FileName>("compression_damage_file"),
           getParam<FileName>("tension_stiffening_file"),
           getParam<FileName>("tension_damage_file"),
           getParam<Real>("youngs_modulus")),
    _local_integrator(_table,
                      {getParam<Real>("youngs_modulus"),
                       getParam<Real>("poissons_ratio"),
                       getParam<Real>("dilation_angle"),
                       getParam<Real>("eccentricity"),
                       getParam<Real>("biaxial_to_uniaxial_compression_ratio"),
                       getParam<Real>("tensile_meridian_ratio"),
                       getParam<unsigned int>("maximum_local_iterations"),
                       getParam<Real>("local_residual_tolerance"),
                       getParam<Real>("local_finite_difference_step"),
                       getParam<Real>("minimum_line_search")}),
    _state_integrator(_local_integrator,
                      {getParam<Real>("tension_recovery"),
                       getParam<Real>("compression_recovery"),
                       getParam<Real>("viscosity"),
                       1.0e-12}),
    _substep_integrator(_state_integrator,
                        {getParam<unsigned int>("maximum_substeps"),
                         getParam<Real>("maximum_strain_increment"),
                         getParam<Real>("reference_tangent_perturbation")}),
    _backbone_plastic_strain(
        declareProperty<RankTwoTensor>(_base_name + "cdp_backbone_plastic_strain")),
    _backbone_plastic_strain_old(
        getMaterialPropertyOld<RankTwoTensor>(_base_name + "cdp_backbone_plastic_strain")),
    _tensile_equivalent_plastic_strain(declareProperty<Real>(_base_name + "cdp_kappa_t")),
    _tensile_equivalent_plastic_strain_old(
        getMaterialPropertyOld<Real>(_base_name + "cdp_kappa_t")),
    _compressive_equivalent_plastic_strain(declareProperty<Real>(_base_name + "cdp_kappa_c")),
    _compressive_equivalent_plastic_strain_old(
        getMaterialPropertyOld<Real>(_base_name + "cdp_kappa_c")),
    _viscous_plastic_strain(
        declareProperty<RankTwoTensor>(_base_name + "cdp_viscous_plastic_strain")),
    _viscous_plastic_strain_old(
        getMaterialPropertyOld<RankTwoTensor>(_base_name + "cdp_viscous_plastic_strain")),
    _damage_t(declareProperty<Real>(_base_name + "DamageT")),
    _damage_t_old(getMaterialPropertyOld<Real>(_base_name + "DamageT")),
    _damage_c(declareProperty<Real>(_base_name + "DamageC")),
    _damage_c_old(getMaterialPropertyOld<Real>(_base_name + "DamageC")),
    _combined_damage(declareProperty<Real>(_base_name + "cdp_combined_damage")),
    _stiffness_factor(declareProperty<Real>(_base_name + "cdp_stiffness_factor")),
    _local_iterations(declareProperty<Real>(_base_name + "cdp_local_iterations")),
    _accepted_substeps(declareProperty<Real>(_base_name + "cdp_accepted_substeps"))
{
}

void
AbaqusCDPStressUpdate::initQpStatefulProperties()
{
  _backbone_plastic_strain[_qp].zero();
  _tensile_equivalent_plastic_strain[_qp] = 0.0;
  _compressive_equivalent_plastic_strain[_qp] = 0.0;
  _viscous_plastic_strain[_qp].zero();
  _damage_t[_qp] = 0.0;
  _damage_c[_qp] = 0.0;
  _combined_damage[_qp] = 0.0;
  _stiffness_factor[_qp] = 1.0;
  _local_iterations[_qp] = 0.0;
  _accepted_substeps[_qp] = 1.0;
}

void
AbaqusCDPStressUpdate::propagateQpStatefulProperties()
{
  _backbone_plastic_strain[_qp] = _backbone_plastic_strain_old[_qp];
  _tensile_equivalent_plastic_strain[_qp] = _tensile_equivalent_plastic_strain_old[_qp];
  _compressive_equivalent_plastic_strain[_qp] = _compressive_equivalent_plastic_strain_old[_qp];
  _viscous_plastic_strain[_qp] = _viscous_plastic_strain_old[_qp];
  _damage_t[_qp] = _damage_t_old[_qp];
  _damage_c[_qp] = _damage_c_old[_qp];
}

TangentCalculationMethod
AbaqusCDPStressUpdate::getTangentCalculationMethod()
{
  return TangentCalculationMethod::FULL;
}

AbaqusCDPStressUpdate::SymmetricTensor
AbaqusCDPStressUpdate::toArray(const RankTwoTensor & tensor)
{
  return {tensor(0, 0), tensor(1, 1), tensor(2, 2), tensor(0, 1), tensor(1, 2), tensor(0, 2)};
}

RankTwoTensor
AbaqusCDPStressUpdate::toRankTwo(const SymmetricTensor & tensor)
{
  RankTwoTensor result;
  result.zero();
  result(0, 0) = tensor[0];
  result(1, 1) = tensor[1];
  result(2, 2) = tensor[2];
  result(0, 1) = result(1, 0) = tensor[3];
  result(1, 2) = result(2, 1) = tensor[4];
  result(0, 2) = result(2, 0) = tensor[5];
  return result;
}

void
AbaqusCDPStressUpdate::assignTangent(const AbaqusCDPSubstepIntegrator::TangentMatrix & source,
                                     RankFourTensor & destination)
{
  destination.zero();
  const std::array<std::array<unsigned int, 2>, 6> pairs = {
      {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2}}};
  for (std::size_t column = 0; column < 6; ++column)
    for (std::size_t row = 0; row < 6; ++row)
    {
      const auto i = pairs[row][0];
      const auto j = pairs[row][1];
      const auto k = pairs[column][0];
      const auto l = pairs[column][1];
      const double value = source[column][row] * (k == l ? 1.0 : 0.5);
      destination(i, j, k, l) = value;
      destination(j, i, k, l) = value;
      destination(i, j, l, k) = value;
      destination(j, i, l, k) = value;
    }
}

AbaqusCDPStressUpdate::CoreState
AbaqusCDPStressUpdate::oldState() const
{
  CoreState result;
  result.backbone.plastic_strain = toArray(_backbone_plastic_strain_old[_qp]);
  result.backbone.tensile_equivalent_plastic_strain = _tensile_equivalent_plastic_strain_old[_qp];
  result.backbone.compressive_equivalent_plastic_strain =
      _compressive_equivalent_plastic_strain_old[_qp];
  result.viscous_plastic_strain = toArray(_viscous_plastic_strain_old[_qp]);
  result.viscous_tension_damage = _damage_t_old[_qp];
  result.viscous_compression_damage = _damage_c_old[_qp];
  return result;
}

void
AbaqusCDPStressUpdate::storeState(const AbaqusCDPSubstepIntegrator::LinearizedResult & result)
{
  const auto & state = result.result.final_result.state;
  _backbone_plastic_strain[_qp] = toRankTwo(state.backbone.plastic_strain);
  _tensile_equivalent_plastic_strain[_qp] = state.backbone.tensile_equivalent_plastic_strain;
  _compressive_equivalent_plastic_strain[_qp] =
      state.backbone.compressive_equivalent_plastic_strain;
  _viscous_plastic_strain[_qp] = toRankTwo(state.viscous_plastic_strain);
  _damage_t[_qp] = state.viscous_tension_damage;
  _damage_c[_qp] = state.viscous_compression_damage;
  _combined_damage[_qp] = result.result.final_result.damage.damage;
  _stiffness_factor[_qp] = result.result.final_result.damage.stiffness_factor;
  _local_iterations[_qp] = result.result.total_local_iterations;
  _accepted_substeps[_qp] = result.result.accepted_substeps;
}

void
AbaqusCDPStressUpdate::updateState(RankTwoTensor & strain_increment,
                                   RankTwoTensor & inelastic_strain_increment,
                                   const RankTwoTensor & /*rotation_increment*/,
                                   RankTwoTensor & stress_new,
                                   const RankTwoTensor & /*stress_old*/,
                                   const RankFourTensor & /*elasticity_tensor*/,
                                   const RankTwoTensor & elastic_strain_old,
                                   const bool compute_full_tangent_operator,
                                   RankFourTensor & tangent_operator)
{
  const CoreState old_state = oldState();
  const auto old_elastic_strain = toArray(elastic_strain_old);
  SymmetricTensor old_total_strain;
  SymmetricTensor new_total_strain;
  const auto increment = toArray(strain_increment);
  for (std::size_t i = 0; i < 6; ++i)
  {
    old_total_strain[i] = old_elastic_strain[i] + old_state.viscous_plastic_strain[i];
    new_total_strain[i] = old_total_strain[i] + increment[i];
  }

  try
  {
    const auto result =
        _substep_integrator.integrateLinearized(old_total_strain, new_total_strain, _dt, old_state);
    stress_new = toRankTwo(result.result.final_result.cauchy_stress);
    SymmetricTensor inelastic_increment_array;
    for (std::size_t i = 0; i < 6; ++i)
      inelastic_increment_array[i] = result.result.final_result.state.viscous_plastic_strain[i] -
                                     old_state.viscous_plastic_strain[i];
    inelastic_strain_increment = toRankTwo(inelastic_increment_array);
    strain_increment -= inelastic_strain_increment;
    if (compute_full_tangent_operator)
      assignTangent(result.algorithmic_tangent, tangent_operator);
    storeState(result);
  }
  catch (const std::exception & error)
  {
    throw MooseException(error.what());
  }
}
