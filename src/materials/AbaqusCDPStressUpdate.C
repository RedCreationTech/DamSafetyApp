#include "AbaqusCDPStressUpdate.h"

#include "MooseException.h"

#include <array>
#include <chrono>
#include <exception>
#include <algorithm>
#include <cctype>
#include <limits>
#include <sstream>

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
  params.addParam<bool>(
      "defer_viscous_update_to_global_step",
      false,
      "Diagnostic-only state-integration switch. Material substeps still advance the inviscid "
      "plasticity backbone and kappa histories, but Duvaut-Lions plastic strain and damage are "
      "committed once using the full accepted global time increment. The default preserves the "
      "existing per-material-substep viscous update.");
  params.addParam<bool>("enable_performance_diagnostics",
                        false,
                        "Measure per-material-call elapsed time and expose detailed local solver "
                        "counters as material properties");
  params.addParam<bool>("enable_path_diagnostics",false,"Opt-in complete costs and bounded path/tangent diagnostics");
  params.addParam<bool>(
      "use_damage_t_secant_tangent",
      false,
      "Diagnostic-only tangent switch. When enabled, a material evaluation that increases "
      "DamageT returns cdp_stiffness_factor times the undamaged elastic tensor to the global "
      "Newton solve. The stress and all constitutive state updates remain unchanged.");
  params.addRangeCheckedParam<Real>(
      "maximum_tensile_history_increment",
      0.0,
      "maximum_tensile_history_increment >= 0",
      "Opt-in upper bound for the accepted cdp_kappa_t increment. A positive value exposes a "
      "material_timestep_limit for IterationAdaptiveDT/MaterialTimeStepPostprocessor without "
      "changing the constitutive update.");
  params.addRangeCheckedParam<Real>(
      "minimum_state_timestep_limit",
      1.0e-8,
      "minimum_state_timestep_limit > 0",
      "Lower bound for the opt-in tensile-history timestep estimate");
  params.addParam<std::vector<unsigned int>>("diagnostic_trace_elements",{},"MOOSE zero-based element IDs for bounded substep traces");
  params.addParam<Real>("diagnostic_time_begin",0.015,"First time for selected path/tangent samples");
  params.addParam<Real>("diagnostic_time_end",0.05,"Last time for selected path/tangent samples");
  params.addParam<unsigned int>("diagnostic_max_trace_calls",2,"Maximum traced material calls per rank/thread object");
  params.addParam<unsigned int>("diagnostic_max_tangent_checks",2,"Maximum plastic tangent samples per rank/thread object");
  params.addParam<unsigned int>("diagnostic_max_failure_samples",4,"Maximum local-failure Jacobian samples per rank/thread object");
  return params;
}

AbaqusCDPStressUpdate::AbaqusCDPStressUpdate(const InputParameters & parameters)
  : StressUpdateBase(parameters),
    _enable_performance_diagnostics(getParam<bool>("enable_performance_diagnostics")),
    _enable_path_diagnostics(getParam<bool>("enable_path_diagnostics")),
    _use_damage_t_secant_tangent(getParam<bool>("use_damage_t_secant_tangent")),
    _maximum_tensile_history_increment(
        getParam<Real>("maximum_tensile_history_increment")),
    _minimum_state_timestep_limit(getParam<Real>("minimum_state_timestep_limit")),
    _state_timestep_limit(std::numeric_limits<Real>::max()),
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
                         getParam<Real>("reference_tangent_perturbation"),
                         getParam<bool>("defer_viscous_update_to_global_step")}),
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
    _damage_t_secant_tangent_active(
        declareProperty<Real>(_base_name + "cdp_damage_t_secant_tangent_active")),
    _local_iterations(declareProperty<Real>(_base_name + "cdp_local_iterations")),
    _jacobian_fallbacks(declareProperty<Real>(_base_name + "cdp_jacobian_fallbacks")),
    _accepted_substeps(declareProperty<Real>(_base_name + "cdp_accepted_substeps")),
    _failed_material_calls(declareProperty<Real>(_base_name + "cdp_failed_material_calls")),
    _attempted_partitions(declareProperty<Real>(_base_name + "cdp_attempted_partitions")),
    _maximum_partition_depth(
        declareProperty<Real>(_base_name + "cdp_maximum_partition_depth")),
    _automatic_jacobian_evaluations(
        declareProperty<Real>(_base_name + "cdp_automatic_jacobian_evaluations")),
    _finite_difference_jacobian_evaluations(
        declareProperty<Real>(_base_name + "cdp_finite_difference_jacobian_evaluations")),
    _local_factorizations(declareProperty<Real>(_base_name + "cdp_local_factorizations")),
    _local_backsolves(declareProperty<Real>(_base_name + "cdp_local_backsolves")),
    _integration_microseconds(
        declareProperty<Real>(_base_name + "cdp_integration_microseconds"))
{
  if (_enable_path_diagnostics)
  {
    std::string label=name();
    for(char & c:label) if(!std::isalnum(static_cast<unsigned char>(c)) && c!='_') c='_';
    const auto suffix=label+"_rank"+std::to_string(processor_id())+"_thread"+std::to_string(_tid);
    _diagnostic_cost_file="cdp_cost_"+suffix+".csv";
    _diagnostic_trace.open("logs/cdp_trace_"+suffix+".jsonl");
    if(!_diagnostic_trace) mooseError("Cannot open CDP diagnostic trace under Job logs/");
  }
}

AbaqusCDPStressUpdate::~AbaqusCDPStressUpdate()
{
  if (_enable_path_diagnostics) writeCostSummary();
}

void
AbaqusCDPStressUpdate::writeCostSummary()
{
  std::ofstream out(_diagnostic_cost_file);
  out<<"category,calls,failed_calls,inclusive_microseconds,failed_inclusive_microseconds\n";
  for(std::size_t i=0;i<CDPDiagnostics::COUNT;++i)
  {
    const auto & c=_diagnostic_cost[i];
    out<<CDPDiagnostics::names[i]<<','<<c.calls<<','<<c.failed<<','<<std::setprecision(17)
       <<c.microseconds<<','<<c.failed_microseconds<<'\n';
  }
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
  _damage_t_secant_tangent_active[_qp] = 0.0;
  _local_iterations[_qp] = 0.0;
  _jacobian_fallbacks[_qp] = 0.0;
  _accepted_substeps[_qp] = 1.0;
  _failed_material_calls[_qp] = 0.0;
  _attempted_partitions[_qp] = 1.0;
  _maximum_partition_depth[_qp] = 0.0;
  _automatic_jacobian_evaluations[_qp] = 0.0;
  _finite_difference_jacobian_evaluations[_qp] = 0.0;
  _local_factorizations[_qp] = 0.0;
  _local_backsolves[_qp] = 0.0;
  _integration_microseconds[_qp] = 0.0;
  _state_timestep_limit = std::numeric_limits<Real>::max();
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
  _damage_t_secant_tangent_active[_qp] = 0.0;
  _state_timestep_limit = std::numeric_limits<Real>::max();
}

Real
AbaqusCDPStressUpdate::computeTimeStepLimit()
{
  return _state_timestep_limit;
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
AbaqusCDPStressUpdate::storeState(const AbaqusCDPSubstepIntegrator::LinearizedResult & result,
                                  const Real integration_microseconds)
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
  _jacobian_fallbacks[_qp] = result.result.total_jacobian_fallbacks;
  _accepted_substeps[_qp] = result.result.accepted_substeps;
  _failed_material_calls[_qp] = result.result.cutback_count;
  _attempted_partitions[_qp] = result.result.attempted_partitions;
  unsigned int partition_depth = 0;
  for (unsigned int substeps = result.result.accepted_substeps; substeps > 1; substeps /= 2)
    ++partition_depth;
  _maximum_partition_depth[_qp] = partition_depth;
  _automatic_jacobian_evaluations[_qp] =
      result.result.total_automatic_jacobian_evaluations;
  _finite_difference_jacobian_evaluations[_qp] =
      result.result.total_finite_difference_jacobian_evaluations;
  _local_factorizations[_qp] = result.result.total_local_factorizations;
  _local_backsolves[_qp] = result.result.total_local_backsolves;
  _integration_microseconds[_qp] = integration_microseconds;
  _state_timestep_limit = std::numeric_limits<Real>::max();
  if (_maximum_tensile_history_increment > 0.0 && _dt > 0.0)
  {
    const Real tensile_history_increment = std::max(
        0.0,
        state.backbone.tensile_equivalent_plastic_strain -
            _tensile_equivalent_plastic_strain_old[_qp]);
    if (tensile_history_increment > _maximum_tensile_history_increment)
      _state_timestep_limit =
          std::max(_minimum_state_timestep_limit,
                   _dt * _maximum_tensile_history_increment / tensile_history_increment);
  }
}

void
AbaqusCDPStressUpdate::updateState(RankTwoTensor & strain_increment,
                                   RankTwoTensor & inelastic_strain_increment,
                                   const RankTwoTensor & /*rotation_increment*/,
                                   RankTwoTensor & stress_new,
                                   const RankTwoTensor & /*stress_old*/,
                                   const RankFourTensor & elasticity_tensor,
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

  CDPDiagnostics::Context context;
  context.counters=&_diagnostic_cost;context.stream=&_diagnostic_trace;
  context.failure_samples=&_failure_samples;
  context.maximum_failure_samples=getParam<unsigned int>("diagnostic_max_failure_samples");
  context.time=_t;context.dt=_dt;context.step=_t_step;
  context.element=_current_elem->id();context.qp=_qp;context.call=++_diagnostic_calls;
  const auto & elements=getParam<std::vector<unsigned int>>("diagnostic_trace_elements");
  const bool selected=_enable_path_diagnostics && _t>=getParam<Real>("diagnostic_time_begin") &&
      _t<=getParam<Real>("diagnostic_time_end") &&
      std::find(elements.begin(),elements.end(),context.element)!=elements.end();
  context.trace=selected && _trace_calls<getParam<unsigned int>("diagnostic_max_trace_calls");
  if(context.trace) ++_trace_calls;
  CDPDiagnostics::Binding binding(_enable_path_diagnostics ? &context : nullptr);
  if(context.trace)
  {
    std::ostringstream payload;
    payload<<"\"old_strain\":";CDPDiagnostics::json(payload,old_total_strain);
    payload<<",\"new_strain\":";CDPDiagnostics::json(payload,new_total_strain);
    payload<<",\"old_state\":";CDPDiagnostics::stateJson(payload,old_state);
    CDPDiagnostics::event("material_input",payload.str());
  }
  try
  {
    const auto integration_start = std::chrono::steady_clock::now();
    const auto result = [&]() {
      CDPDiagnostics::Scope material_scope(CDPDiagnostics::MATERIAL);
      return _substep_integrator.integrateLinearized(old_total_strain,new_total_strain,_dt,old_state);
    }();
    const Real integration_microseconds =
        _enable_performance_diagnostics
            ? std::chrono::duration<Real, std::micro>(std::chrono::steady_clock::now() -
                                                       integration_start)
                  .count()
            : 0.0;
    if(selected && result.result.final_result.backbone.plastic &&
        _tangent_checks<getParam<unsigned int>("diagnostic_max_tangent_checks"))
    {
      ++_tangent_checks;
      auditTangent(old_total_strain,new_total_strain,old_state,result);
    }
    if(_enable_path_diagnostics && _diagnostic_calls%1000==0) writeCostSummary();
    stress_new = toRankTwo(result.result.final_result.cauchy_stress);
    SymmetricTensor inelastic_increment_array;
    for (std::size_t i = 0; i < 6; ++i)
      inelastic_increment_array[i] = result.result.final_result.state.viscous_plastic_strain[i] -
                                     old_state.viscous_plastic_strain[i];
    inelastic_strain_increment = toRankTwo(inelastic_increment_array);
    strain_increment -= inelastic_strain_increment;
    const bool damage_t_secant_tangent_active =
        _use_damage_t_secant_tangent &&
        result.result.final_result.state.viscous_tension_damage >
            old_state.viscous_tension_damage + 1.0e-12;
    _damage_t_secant_tangent_active[_qp] = damage_t_secant_tangent_active ? 1.0 : 0.0;
    if (compute_full_tangent_operator)
    {
      if (damage_t_secant_tangent_active)
      {
        tangent_operator = elasticity_tensor;
        tangent_operator *= result.result.final_result.damage.stiffness_factor;
      }
      else
        assignTangent(result.algorithmic_tangent, tangent_operator);
    }
    storeState(result, integration_microseconds);
  }
  catch (const std::exception & error)
  {
    if(_enable_path_diagnostics)
    {
      std::ostringstream payload;
      payload<<"\"old_strain\":";CDPDiagnostics::json(payload,old_total_strain);
      payload<<",\"new_strain\":";CDPDiagnostics::json(payload,new_total_strain);
      payload<<",\"old_state\":";CDPDiagnostics::stateJson(payload,old_state);
      CDPDiagnostics::event("material_failure",payload.str());
      writeCostSummary();
    }
    throw MooseException(error.what());
  }
}

void
AbaqusCDPStressUpdate::auditTangent(
    const SymmetricTensor & old_strain,const SymmetricTensor & new_strain,
    const CoreState & old_state,const AbaqusCDPSubstepIntegrator::LinearizedResult & result)
{
  std::ostringstream payload;
  payload<<"\"old_strain\":";CDPDiagnostics::json(payload,old_strain);
  payload<<",\"new_strain\":";CDPDiagnostics::json(payload,new_strain);
  payload<<",\"old_state\":";CDPDiagnostics::stateJson(payload,old_state);
  payload<<",\"base_partition\":"<<result.result.accepted_substeps;
  payload<<",\"algorithmic_tangent_columns\":";CDPDiagnostics::json(payload,result.algorithmic_tangent);
  payload<<",\"samples\":[";
  bool first=true;
  // No diagnostic recomputation is counted as a production material call.
  {
    CDPDiagnostics::Binding exclude_diagnostic_cost(nullptr);
    for(double h : {1e-8,1e-9,1e-10})
      for(std::size_t column=0;column<6;++column)
      {
        if(!first)payload<<',';first=false;
        payload<<"{\"h\":"<<h<<",\"column\":"<<column;
        auto plus=new_strain,minus=new_strain;plus[column]+=h;minus[column]-=h;
        try
        {
          const auto rp=_substep_integrator.integrateLinearized(old_strain,plus,_dt,old_state);
          const auto rm=_substep_integrator.integrateLinearized(old_strain,minus,_dt,old_state);
          double diff2=0,norm2=0;
          SymmetricTensor centered={},forward={},backward={};
          for(std::size_t row=0;row<6;++row)
          {
            const double p=rp.result.final_result.cauchy_stress[row];
            const double m=rm.result.final_result.cauchy_stress[row];
            const double b=result.result.final_result.cauchy_stress[row];
            centered[row]=(p-m)/(2*h);forward[row]=(p-b)/h;backward[row]=(b-m)/h;
            const double error=centered[row]-result.algorithmic_tangent[column][row];
            diff2+=error*error;norm2+=centered[row]*centered[row];
          }
          payload<<",\"ok\":true,\"plus_partition\":"<<rp.result.accepted_substeps
                 <<",\"minus_partition\":"<<rm.result.accepted_substeps
                 <<",\"same_final_branch\":"
                 <<((rp.result.final_result.backbone.active_branch==result.result.final_result.backbone.active_branch &&
                     rm.result.final_result.backbone.active_branch==result.result.final_result.backbone.active_branch)?"true":"false")
                 <<",\"relative_error\":";
          CDPDiagnostics::json(payload,std::sqrt(diff2)/std::max(std::sqrt(norm2),getParam<Real>("youngs_modulus")*1e-12));
          payload<<",\"centered\":";CDPDiagnostics::json(payload,centered);
          payload<<",\"forward\":";CDPDiagnostics::json(payload,forward);
          payload<<",\"backward\":";CDPDiagnostics::json(payload,backward);
        }
        catch(const std::exception &) { payload<<",\"ok\":false"; }
        payload<<'}';
      }
    // Repeat exactly the original call after perturbations: old state must be unchanged.
    try
    {
      const auto repeat=_substep_integrator.integrateLinearized(old_strain,new_strain,_dt,old_state);
      double error=0;
      for(std::size_t i=0;i<6;++i)
        error=std::max(error,std::abs(repeat.result.final_result.cauchy_stress[i]-result.result.final_result.cauchy_stress[i]));
      payload<<"],\"repeat_max_stress_error\":";CDPDiagnostics::json(payload,error);
      payload<<",\"repeat_same_partition\":"<<(repeat.result.accepted_substeps==result.result.accepted_substeps?"true":"false");
    }
    catch(const std::exception &) { payload<<"],\"repeat_failed\":true"; }
  }
  CDPDiagnostics::event("material_tangent",payload.str());
}
