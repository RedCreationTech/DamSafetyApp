#include "CDPAssemblyProbe.h"
#include "libmesh/elem.h"
#include <fstream>
#include <iomanip>

registerMooseObject("DamSafetyApp", CDPAssemblyProbe);
namespace
{
std::ofstream cdp_probe_stream;
bool cdp_probe_tangent = false;
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
  for (const auto id : getParam<std::vector<dof_id_type>>("probe_elements"))
    _probe_elements.insert(id);
}
void CDPAssemblyProbe::beginCapture(const std::string & path, bool tangent)
{
  if (libMesh::n_threads() != 1)
    mooseError("CDPAssemblyProbe requires one thread per MPI rank");
  cdp_probe_tangent = tangent;
  cdp_probe_stream.open(path);
  if (!cdp_probe_stream) mooseError("Cannot open assembly probe capture");
  cdp_probe_stream << std::setprecision(17) << "element,qp,x,y,z,kappa_t,DamageT,substeps,depth,fallbacks";
  for (const auto name : {"stress", "strain"})
    for (unsigned int c=0;c<9;++c) cdp_probe_stream << ',' << name << c;
  if (tangent)
    for (unsigned int c=0;c<81;++c) cdp_probe_stream << ",C" << c;
  cdp_probe_stream << '\n';
}
void CDPAssemblyProbe::endCapture()
{
  if (!cdp_probe_stream) mooseError("Assembly probe write failed");
  cdp_probe_stream.close();
}
void CDPAssemblyProbe::computeQpProperties()
{
  if (!cdp_probe_stream.is_open() || !_probe_elements.count(_current_elem->id())) return;
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
  cdp_probe_stream << '\n';
}
