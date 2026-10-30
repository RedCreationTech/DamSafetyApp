#include "AqueductElementSampler.h"
#include "MooseMesh.h"
#include "libmesh/elem.h"

registerMooseObject("DamSafetyApp",AqueductElementSampler);
InputParameters AqueductElementSampler::validParams()
{
  auto p=ElementValueSampler::validParams();
  p.addRequiredParam<std::vector<dof_id_type>>("element_ids","Native zero-based source-mapped element IDs");
  return p;
}
AqueductElementSampler::AqueductElementSampler(const InputParameters & p)
 :ElementValueSampler(p),_ids(getParam<std::vector<dof_id_type>>("element_ids").begin(),
                             getParam<std::vector<dof_id_type>>("element_ids").end())
{
  if(_ids.empty() || _ids.size()!=getParam<std::vector<dof_id_type>>("element_ids").size())
    paramError("element_ids","A nonempty list of unique element IDs is required");
}
void AqueductElementSampler::initialSetup()
{
  ElementValueSampler::initialSetup();
  // The full aqueduct uses a replicated mesh. Refuse missing IDs or block omissions.
  for(auto id:_ids)
  {
    const auto * e=_mesh.getMesh().query_elem_ptr(id);
    if(!e || !e->active() || !hasBlocks(e->subdomain_id()))
      paramError("element_ids","Element missing, inactive or outside sampler blocks: ",id);
  }
}
void AqueductElementSampler::execute()
{
  if(_ids.count(_current_elem->id()))ElementValueSampler::execute();
}
