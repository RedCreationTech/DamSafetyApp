#include "CDPAssemblyProbe.h"
#include "libmesh/elem.h"
#include <fstream>
#include <iomanip>

registerMooseObject("DamSafetyApp", CDPAssemblyProbe);
namespace
{
std::ofstream cdp_probe_stream;
std::ofstream cdp_substep_stream;
unsigned int cdp_substep_call = 0;
bool cdp_probe_tangent = false;
bool cdp_probe_all = false;
const std::vector<std::string> tensor_history_names = {
    "cdp_backbone_plastic_strain", "cdp_viscous_plastic_strain",
    "elastic_strain", "mechanical_strain"};
const std::vector<std::string> scalar_history_names = {
    "cdp_kappa_t", "cdp_kappa_c", "DamageT", "DamageC"};
}
InputParameters CDPAssemblyProbe::validParams()
{
  auto p = Material::validParams();
  p.addRequiredParam<std::vector<dof_id_type>>("probe_elements", "Selected libMesh element IDs");
  p.addClassDescription("Passive material observer gated by accepted-state Output");
  return p;
}
CDPAssemblyProbe::CDPAssemblyProbe(const InputParameters & p)
  : Material(p),
    _probe_tangent(getMaterialProperty<RankFourTensor>("Jacobian_mult")),
    _probe_stress(getMaterialProperty<RankTwoTensor>("stress")),
    _probe_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")),
    _probe_kappa(getMaterialProperty<Real>("cdp_kappa_t")),
    _probe_damage(getMaterialProperty<Real>("DamageT")),
    _probe_substeps(getMaterialProperty<Real>("cdp_accepted_substeps")),
    _probe_depth(getMaterialProperty<Real>("cdp_maximum_partition_depth")),
    _probe_fallbacks(getMaterialProperty<Real>("cdp_jacobian_fallbacks"))
{
  for (const auto & name : tensor_history_names)
  {
    _history_tensors.push_back(&getMaterialProperty<RankTwoTensor>(name));
    _history_tensors.push_back(&getMaterialPropertyOld<RankTwoTensor>(name));
  }
  for (const auto & name : scalar_history_names)
  {
    _history_scalars.push_back(&getMaterialProperty<Real>(name));
    _history_scalars.push_back(&getMaterialPropertyOld<Real>(name));
  }
  _history_scalars.push_back(&getMaterialProperty<Real>("cdp_stiffness_factor"));
  for (const auto & name : {"cdp_audit_branch", "cdp_audit_residual", "cdp_audit_plastic"})
    _history_scalars.push_back(&getMaterialProperty<Real>(name));
  for (const auto id : getParam<std::vector<dof_id_type>>("probe_elements"))
    _probe_elements.insert(id);
}
void CDPAssemblyProbe::beginCapture(const std::string & path, bool tangent, bool all_elements, bool substeps)
{
  if (libMesh::n_threads() != 1)
    ::mooseError("CDPAssemblyProbe requires one thread per MPI rank");
  if (substeps)
  {
    cdp_substep_call = 0;
    cdp_substep_stream.open(path + "_substeps.csv");
    if (!cdp_substep_stream) ::mooseError("Cannot open passive CDP substep trace");
    cdp_substep_stream << std::setprecision(17)
                       << "element,qp,call,partition,substep,substep_dt,succeeded,error";
    for (const auto name : {"target", "stress", "effective_stress"})
      for (unsigned int i = 0; i < 6; ++i) cdp_substep_stream << ',' << name << i;
    for (const auto age : {"old", "new"})
      for (unsigned int i = 0; i < 16; ++i) cdp_substep_stream << ',' << age << "_state" << i;
    cdp_substep_stream << ",branch,plastic,iterations,residual_norm,trial_yield,final_yield,plastic_multiplier"
                       << ",jacobian_fallbacks,automatic_jacobian_evaluations,finite_difference_jacobian_evaluations";
    for (unsigned int i = 0; i < 22 * 22; ++i) cdp_substep_stream << ",transition" << i;
    for (unsigned int i = 0; i < 36; ++i) cdp_substep_stream << ",chained" << i;
    cdp_substep_stream << '\n';
  }
  cdp_probe_tangent = tangent;
  cdp_probe_all = all_elements;
  cdp_probe_stream.open(path);
  if (!cdp_probe_stream) ::mooseError("Cannot open assembly probe capture");
  cdp_probe_stream << std::setprecision(17) << "element,qp,x,y,z,kappa_t,DamageT,substeps,depth,fallbacks";
  for (const auto name : {"stress", "strain"})
    for (unsigned int c=0;c<9;++c) cdp_probe_stream << ',' << name << c;
  if (tangent)
    for (unsigned int c=0;c<81;++c) cdp_probe_stream << ",C" << c;
  for (const auto & name : tensor_history_names)
    for (const auto age : {"new", "old"})
      for (unsigned int c = 0; c < 9; ++c)
        cdp_probe_stream << ',' << name << '_' << age << c;
  for (const auto & name : scalar_history_names)
    for (const auto age : {"new", "old"})
      cdp_probe_stream << ',' << name << '_' << age;
  cdp_probe_stream << ",cdp_stiffness_factor,final_substep_active_branch,final_substep_residual_norm,final_substep_plastic\n";
}
void CDPAssemblyProbe::endCapture()
{
  if (!cdp_probe_stream) ::mooseError("Assembly probe write failed");
  cdp_probe_stream.close();
  if (cdp_substep_stream.is_open())
  {
    if (!cdp_substep_stream) ::mooseError("Passive CDP substep trace write failed");
    cdp_substep_stream.close();
  }
}
void CDPAssemblyProbe::computeQpProperties()
{
  if (!cdp_probe_stream.is_open() || (!cdp_probe_all && !_probe_elements.count(_current_elem->id()))) return;
  cdp_probe_stream << _current_elem->id() << ',' << _qp;
  for (unsigned int c=0;c<3;++c) cdp_probe_stream << ',' << _q_point[_qp](c);
  cdp_probe_stream << ',' << _probe_kappa[_qp] << ',' << _probe_damage[_qp]
                   << ',' << _probe_substeps[_qp] << ',' << _probe_depth[_qp]
                   << ',' << _probe_fallbacks[_qp];
  for (const auto * tensor : {&_probe_stress[_qp], &_probe_strain[_qp]})
    for (unsigned int c=0;c<9;++c) cdp_probe_stream << ',' << (*tensor)(c/3,c%3);
  if (cdp_probe_tangent)
    for (unsigned int c=0;c<81;++c)
      cdp_probe_stream << ',' << _probe_tangent[_qp](c/27,(c/9)%3,(c/3)%3,c%3);
  for (const auto * property : _history_tensors)
    for (unsigned int c = 0; c < 9; ++c)
      cdp_probe_stream << ',' << (*property)[_qp](c / 3, c % 3);
  for (const auto * property : _history_scalars)
    cdp_probe_stream << ',' << (*property)[_qp];
  cdp_probe_stream << '\n';
}

bool CDPAssemblyProbe::capturingSubsteps()
{
  return cdp_substep_stream.is_open();
}
void CDPAssemblyProbe::recordSubsteps(dof_id_type element, unsigned int qp,
                                     const AbaqusCDPSubstepIntegrator::Trace & trace)
{
  if (!capturingSubsteps()) return;
  const auto call = ++cdp_substep_call;
  auto state = [](const AbaqusCDPSubstepIntegrator::State & s) {
    for (const auto x : s.backbone.plastic_strain) cdp_substep_stream << ',' << x;
    cdp_substep_stream << ',' << s.backbone.tensile_equivalent_plastic_strain
                       << ',' << s.backbone.compressive_equivalent_plastic_strain;
    for (const auto x : s.viscous_plastic_strain) cdp_substep_stream << ',' << x;
    cdp_substep_stream << ',' << s.viscous_tension_damage << ',' << s.viscous_compression_damage;
  };
  for (const auto & e : trace)
  {
    cdp_substep_stream << element << ',' << qp << ',' << call << ',' << e.partition << ','
                       << e.substep << ',' << e.substep_dt << ',' << e.succeeded << ",\"";
    for (const char c : e.error)
    {
      if (c == '"') cdp_substep_stream << '"';
      cdp_substep_stream << c;
    }
    cdp_substep_stream << '"';
    for (const auto x : e.target) cdp_substep_stream << ',' << x;
    // A failed step has no converged output/derivative; leave these cells empty.
    for (const auto * t : {&e.result.cauchy_stress, &e.result.viscous_effective_stress})
      for (const auto x : *t) { cdp_substep_stream << ','; if (e.succeeded) cdp_substep_stream << x; }
    state(e.old_state);
    if (e.succeeded) state(e.result.state);
    else for (unsigned int i = 0; i < 16; ++i) cdp_substep_stream << ',';
    if (e.succeeded)
    {
      const auto & b = e.result.backbone;
      cdp_substep_stream << ',' << static_cast<unsigned int>(b.active_branch) << ',' << b.plastic
                         << ',' << b.iterations << ',' << b.residual_norm << ',' << b.trial_yield
                         << ',' << b.final_yield << ',' << b.plastic_multiplier << ',' << b.jacobian_fallbacks
                         << ',' << b.automatic_jacobian_evaluations << ',' << b.finite_difference_jacobian_evaluations;
    }
    else for (unsigned int i = 0; i < 10; ++i) cdp_substep_stream << ',';
    for (const auto & column : e.transition)
      for (const auto x : column) { cdp_substep_stream << ','; if (e.succeeded) cdp_substep_stream << x; }
    for (const auto & column : e.chained_tangent)
      for (const auto x : column) { cdp_substep_stream << ','; if (e.succeeded) cdp_substep_stream << x; }
    cdp_substep_stream << '\n';
  }
  if (!cdp_substep_stream) ::mooseError("Passive CDP substep trace write failed");
}
