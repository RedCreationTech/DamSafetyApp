#include "AqueductEmbeddedConstraint.h"
#include "MooseMesh.h"
#include "libmesh/elem.h"
#include <algorithm>
#include <array>
#include <fstream>
#include <sstream>
#include <set>
#include <cmath>

registerMooseObject("DamSafetyApp", AqueductEmbeddedConstraint);

InputParameters AqueductEmbeddedConstraint::validParams()
{
  auto p=EqualValueEmbeddedConstraint::validParams();
  p.addClassDescription("Host-restricted embedded kinematics from an audited unique node-to-host map");
  p.addRequiredParam<FileName>("source_map","ASCII map: node,host,owning_secondary_block,x,y,z,host_node_0..7; runtime IDs");
  p.addRequiredParam<std::vector<SubdomainID>>("allowed_primary_block_ids","Only source-approved host blocks");
  return p;
}
AqueductEmbeddedConstraint::AqueductEmbeddedConstraint(const InputParameters & p)
 :EqualValueEmbeddedConstraint(p)
{
  prepareSecondaryToPrimaryMap();
}
void AqueductEmbeddedConstraint::prepareSecondaryToPrimaryMap()
{
  const auto allowed=getParam<std::vector<SubdomainID>>("allowed_primary_block_ids");
  if(allowed.empty())paramError("allowed_primary_block_ids","Source host blocks required");
  std::ifstream input(getParam<FileName>("source_map"));
  if(!input)paramError("source_map","Cannot read audited embedded map");
  std::string line;std::getline(input,line);if(!line.empty() && line.back()=='\r')line.pop_back();
  if(line!="secondary_node,primary_element,secondary_block,x_m,y_m,z_m,host_node_0,host_node_1,host_node_2,host_node_3,host_node_4,host_node_5,host_node_6,host_node_7")
    paramError("source_map","Invalid map schema");
  _secondary_to_primary_map.clear();std::set<dof_id_type> seen;
  unsigned long row=1;
  while(std::getline(input,line))
  {
    ++row;if(line.empty())continue;std::replace(line.begin(),line.end(),',',' ');
    std::istringstream values(line);dof_id_type sid,eid,owner;Point expected;
    if(!(values>>sid>>eid>>owner>>expected(0)>>expected(1)>>expected(2)))
      paramError("source_map","Malformed row "+std::to_string(row));
    std::array<dof_id_type,8> host_nodes;
    for(auto & node:host_nodes)if(!(values>>node))paramError("source_map","Missing host topology");
    std::string extra;if(values>>extra)paramError("source_map","Unexpected extra fields");
    if(!seen.insert(sid).second)paramError("source_map","Each embedded node must have exactly one owner");
    auto node=_mesh.getMesh().query_node_ptr(sid);
    if(!node)paramError("source_map","Missing secondary node; use the frozen replicated mesh");
    for(unsigned j=0;j<3;++j)
      if(!std::isfinite(expected(j)) || std::abs((*node)(j)-expected(j))>1e-11)
        paramError("source_map","Secondary node identity/coordinate mismatch");
    const auto it=_mesh.nodeToElemMap().find(sid);bool belongs=false;
    if(it!=_mesh.nodeToElemMap().end())
      for(auto child:it->second)
      {
        auto elem=_mesh.elemPtr(child);
        if(elem->subdomain_id()==owner && elem->dim()==1 && elem->n_nodes()==2)belongs=true;
      }
    if(!belongs)paramError("source_map","Owning truss block does not contain the secondary node");
    if(owner!=_secondary)continue;
    auto host=_mesh.getMesh().query_elem_ptr(eid);
    if(!host || host->n_nodes()!=8 || host->dim()!=3)
      paramError("source_map","Source host is not a three-dimensional HEX8");
    if(std::find(allowed.begin(),allowed.end(),host->subdomain_id())==allowed.end())
      paramError("source_map","Host is outside the source-approved blocks");
    for(unsigned j=0;j<8;++j)
      if(host->node_id(j)!=host_nodes[j])paramError("source_map","Host element identity/topology mismatch");
    if(!host->contains_point(*node,1e-8))paramError("source_map","Node lies outside its declared source host");
    _secondary_to_primary_map.emplace(sid,eid);_subproblem.addGhostedElem(eid);
  }
  if(_secondary_to_primary_map.empty())paramError("source_map","This secondary block owns no embedded nodes");
}
