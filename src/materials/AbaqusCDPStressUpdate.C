#include "AbaqusCDPStressUpdate.h"

#include "MooseException.h"
#include "FEProblemBase.h"
#include "libmesh/elem.h"

#include <array>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <limits>

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
  params.addParam<bool>("enable_performance_diagnostics",
                        false,
                        "Measure per-material-call elapsed time and expose detailed local solver "
                        "counters as material properties");
  params.addParam<bool>("enable_exact_result_cache", false,
                        "Opt-in reuse of complete successful results for bit-identical point inputs");
  params.addParam<bool>("observe_exact_result_cache", false,
                        "Count potential cache hits but always execute the original integrator");
  params.addRangeCheckedParam<unsigned int>("exact_result_cache_capacity", 65536,
                                             "exact_result_cache_capacity > 0",
                                             "Maximum cached points per material/thread");
  params.addParam<std::string>("exact_cache_diagnostics_prefix", "",
                               "Optional per-rank/thread lifetime cache counters, written on shutdown");
  params.addParam<std::string>("failure_point_capture_prefix", "",
                               "Empty disables capture; otherwise bounded per-rank/thread point packets");
  params.addRangeCheckedParam<unsigned int>("failure_point_capture_limit", 1,
                                             "failure_point_capture_limit > 0", "Capture attempts per material/thread");
  params.addParam<Real>("failure_point_capture_start", 0., "Total-time window start");
  params.addParam<Real>("failure_point_capture_end", std::numeric_limits<Real>::max(),
                       "Total-time window end");
  return params;
}

AbaqusCDPStressUpdate::AbaqusCDPStressUpdate(const InputParameters & parameters)
  : StressUpdateBase(parameters),
    _enable_performance_diagnostics(getParam<bool>("enable_performance_diagnostics")),
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
    _enable_exact_result_cache(getParam<bool>("enable_exact_result_cache")),
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
    _jacobian_fallbacks(declareProperty<Real>(_base_name + "cdp_jacobian_fallbacks")),
    _audit_branch(declareProperty<Real>(_base_name + "cdp_audit_branch")),
    _audit_residual(declareProperty<Real>(_base_name + "cdp_audit_residual")),
    _audit_plastic(declareProperty<Real>(_base_name + "cdp_audit_plastic")),
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
  if (_enable_exact_result_cache && getParam<bool>("observe_exact_result_cache"))
    paramError("observe_exact_result_cache", "Observe and reuse modes are mutually exclusive");
  if (_enable_exact_result_cache || getParam<bool>("observe_exact_result_cache"))
    _exact_cache = std::make_unique<AbaqusCDPExactResultCache>(
        getParam<unsigned int>("exact_result_cache_capacity"));
  std::string safe_name=name();
  for(char & c:safe_name)if(!std::isalnum(static_cast<unsigned char>(c)))c='_';
  const auto suffix="."+safe_name+".rank"+std::to_string(processor_id())+
                    ".thread"+std::to_string(_tid);
  const auto & diagnostics=getParam<std::string>("exact_cache_diagnostics_prefix");
  if(!diagnostics.empty())
  {
    if(!_exact_cache)paramError("exact_cache_diagnostics_prefix", "Requires observe or reuse mode");
    const auto path=diagnostics+suffix+".csv";
    if(std::filesystem::exists(path))paramError("exact_cache_diagnostics_prefix", "File already exists");
    _cache_diagnostics_path=path;
  }
  const auto start=getParam<Real>("failure_point_capture_start");
  const auto end=getParam<Real>("failure_point_capture_end");
  if(!std::isfinite(start) || !std::isfinite(end) || end<start)
    paramError("failure_point_capture_end", "Capture window must be finite and ordered");
  if(!getParam<std::string>("failure_point_capture_prefix").empty())
  {
    _capture_template=std::make_unique<AbaqusCDPPointPacket>();
    _capture_template->local=_local_integrator.parameters();
    _capture_template->state=_state_integrator.parameters();
    _capture_template->substeps=_substep_integrator.parameters();
    _capture_template->material=name();
    _capture_template->setTables(_table);
  }
}

AbaqusCDPStressUpdate::~AbaqusCDPStressUpdate()
{
  if(_exact_cache && _cache_diagnostics.is_open())
  {
    const auto & s=_exact_cache->counters();
    _cache_diagnostics << s.requests << ',' << s.potential_hits << ',' << s.reused_hits << ','
        << _cache_integrations << ',' << s.stores << ',' << s.capacity_skips << ','
        << s.invalidations << ',' << s.peak_entries << ',' << _exact_cache->entryPayloadBytes() << ','
        << _cache_integrator_microseconds << ',' << _cache_lookup_store_microseconds << '\n';
    _cache_diagnostics.flush();
    if(!_cache_diagnostics)std::fprintf(stderr,"CDP cache diagnostic write failed\n");
  }
}

void AbaqusCDPStressUpdate::timestepSetup()
{
  StressUpdateBase::timestepSetup();
  if(_exact_cache)_exact_cache->invalidate();
}
void AbaqusCDPStressUpdate::residualSetup()
{
  StressUpdateBase::residualSetup();
  // Retain entries only from the latest residual evaluation to its Jacobian.
  if(_exact_cache)_exact_cache->invalidate();
}
void AbaqusCDPStressUpdate::meshChanged()
{
  StressUpdateBase::meshChanged();
  if(_exact_cache)_exact_cache->invalidate();
}

AbaqusCDPSubstepIntegrator::LinearizedResult
AbaqusCDPStressUpdate::integrateWithOptionalCache(const SymmetricTensor & old_strain,
                                                 const SymmetricTensor & new_strain,
                                                 const CoreState & old_state)
{
  if(!_exact_cache)
    return _substep_integrator.integrateLinearized(old_strain,new_strain,_dt,old_state);
  // C06 checks the input before solving in the same directory; checks must not create output.
  if(!_cache_diagnostics_path.empty() && !_cache_diagnostics.is_open())
  {
    if(std::filesystem::exists(_cache_diagnostics_path))
      paramError("exact_cache_diagnostics_prefix", "File already exists");
    _cache_diagnostics.open(_cache_diagnostics_path);
    if(!_cache_diagnostics)paramError("exact_cache_diagnostics_prefix", "Cannot create diagnostics");
    _cache_diagnostics << "requests,potential_hits,reused_hits,integrations,stores,capacity_skips,"
                          "invalidations,peak_entries,entry_payload_bytes,integrator_microseconds,"
                          "lookup_store_microseconds\n";
    _cache_diagnostics.flush();
  }
  const auto lookup_start=std::chrono::steady_clock::now();
  const auto key=AbaqusCDPExactResultCache::key(_current_elem->id(),_qp,old_strain,new_strain,
      _dt,old_state,(_bnd || _neighbor) ? _current_side : std::numeric_limits<unsigned int>::max());
  const auto * hit=_exact_cache->find(key,_enable_exact_result_cache);
  _cache_lookup_store_microseconds += std::chrono::duration<double,std::micro>(
      std::chrono::steady_clock::now()-lookup_start).count();
  if(hit && _enable_exact_result_cache)return *hit;
  ++_cache_integrations;
  const auto core_start=std::chrono::steady_clock::now();
  auto result = [&]() {
    try
    {
      auto computed=_substep_integrator.integrateLinearized(old_strain,new_strain,_dt,old_state);
      _cache_integrator_microseconds += std::chrono::duration<double,std::micro>(
          std::chrono::steady_clock::now()-core_start).count();
      return computed;
    }
    catch(...)
    {
      _cache_integrator_microseconds += std::chrono::duration<double,std::micro>(
          std::chrono::steady_clock::now()-core_start).count();
      throw;
    }
  }();
  const auto store_start=std::chrono::steady_clock::now();
  _exact_cache->store(key,result);
  _cache_lookup_store_microseconds += std::chrono::duration<double,std::micro>(
      std::chrono::steady_clock::now()-store_start).count();
  return result;
}

void AbaqusCDPStressUpdate::captureFailure(const SymmetricTensor & old_strain,
                                           const SymmetricTensor & new_strain,
                                           const CoreState & old_state,
                                           const char * error) noexcept
{
  if(!_capture_template || _t<getParam<Real>("failure_point_capture_start") ||
     _t>getParam<Real>("failure_point_capture_end") ||
     _capture_attempts>=getParam<unsigned int>("failure_point_capture_limit"))return;
  const auto record=_capture_attempts++;
  try
  {
    auto packet=*_capture_template;
    packet.element=_current_elem->id();packet.qp=_qp;packet.rank=processor_id();packet.thread=_tid;
    packet.side=(_bnd || _neighbor) ? _current_side : std::numeric_limits<unsigned int>::max();
    packet.step=_fe_problem.timeStep();packet.time=_t;packet.dt=_dt;
    for(unsigned int j=0;j<3;++j)packet.position[j]=_q_point[_qp](j);
    packet.old_strain=old_strain;packet.new_strain=new_strain;packet.old=old_state;
    packet.expected_error=error;
    std::string safe_name=name();
    for(char & c:safe_name)if(!std::isalnum(static_cast<unsigned char>(c)))c='_';
    const auto path=getParam<std::string>("failure_point_capture_prefix")+"."+safe_name+
        ".rank"+std::to_string(packet.rank)+".thread"+std::to_string(packet.thread)+
        ".record"+std::to_string(record)+".txt";
    if(std::filesystem::exists(path))throw std::runtime_error("Capture file already exists");
    std::ofstream out(path);packet.write(out);out.flush();
    if(!out)throw std::runtime_error("Capture flush failed");
  }
  catch(const std::exception & e)
  {
    // Diagnostic I/O must never replace the original material exception.
    std::fprintf(stderr,"CDP point capture failed: %s\n",e.what());
  }
  catch(...){std::fprintf(stderr,"CDP point capture failed: unknown diagnostic error\n");}
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
  _jacobian_fallbacks[_qp] = 0.0;
  _audit_branch[_qp] = _audit_residual[_qp] = _audit_plastic[_qp] = 0.;
  _accepted_substeps[_qp] = 1.0;
  _failed_material_calls[_qp] = 0.0;
  _attempted_partitions[_qp] = 1.0;
  _maximum_partition_depth[_qp] = 0.0;
  _automatic_jacobian_evaluations[_qp] = 0.0;
  _finite_difference_jacobian_evaluations[_qp] = 0.0;
  _local_factorizations[_qp] = 0.0;
  _local_backsolves[_qp] = 0.0;
  _integration_microseconds[_qp] = 0.0;
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
  _audit_branch[_qp] = static_cast<unsigned int>(result.result.final_result.backbone.active_branch);
  _audit_residual[_qp] = result.result.final_result.backbone.residual_norm;
  _audit_plastic[_qp] = result.result.final_result.backbone.plastic;
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
    const auto integration_start = std::chrono::steady_clock::now();
    const auto result = integrateWithOptionalCache(old_total_strain,new_total_strain,old_state);
    const Real integration_microseconds =
        _enable_performance_diagnostics
            ? std::chrono::duration<Real, std::micro>(std::chrono::steady_clock::now() -
                                                       integration_start)
                  .count()
            : 0.0;
    stress_new = toRankTwo(result.result.final_result.cauchy_stress);
    SymmetricTensor inelastic_increment_array;
    for (std::size_t i = 0; i < 6; ++i)
      inelastic_increment_array[i] = result.result.final_result.state.viscous_plastic_strain[i] -
                                     old_state.viscous_plastic_strain[i];
    inelastic_strain_increment = toRankTwo(inelastic_increment_array);
    strain_increment -= inelastic_strain_increment;
    if (compute_full_tangent_operator)
      assignTangent(result.algorithmic_tangent, tangent_operator);
    storeState(result, integration_microseconds);
  }
  catch (const std::exception & error)
  {
    if(_exact_cache)_exact_cache->invalidate();
    captureFailure(old_total_strain,new_total_strain,old_state,error.what());
    throw MooseException(error.what());
  }
}
