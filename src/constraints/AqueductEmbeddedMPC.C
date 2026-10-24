#include "AqueductEmbeddedMPC.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SystemBase.h"
#include "Assembly.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/sparse_matrix.h"
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <cmath>

registerMooseObject("DamSafetyApp",AqueductEmbeddedMPC);
InputParameters AqueductEmbeddedMPC::validParams()
{
  auto p=NodalConstraint::validParams();
  p.addRequiredCoupledVar("displacements","Three physical translation variables");
  p.addRequiredCoupledVar("multipliers","Three explicit current-time constraint multipliers");
  p.addRequiredParam<FileName>("source_map","Audited unique embedded node/host geometry map");
  p.addRequiredParam<std::vector<SubdomainID>>("allowed_primary_block_ids","Source-approved solid host blocks");
  p.suppressParameter<MooseEnum>("formulation");
  p.suppressParameter<NonlinearVariableName>("variable_secondary");
  p.set<bool>("use_displaced_mesh")=false;
  p.addClassDescription("Batched conservative host-restricted embedded MPC, including dynamic mass transfer");
  return p;
}
AqueductEmbeddedMPC::AqueductEmbeddedMPC(const InputParameters & p):NodalConstraint(p)
{
  if(getParam<bool>("use_displaced_mesh"))paramError("use_displaced_mesh","Source SMALL strain uses reference coordinates");
  if(coupledComponents("displacements")!=3 || coupledComponents("multipliers")!=3)
    paramError("multipliers","Exactly three physical variables and three multipliers required");
  for(unsigned j=0;j<3;++j)
  {
    _physical[j]=getVar("displacements",j);_multipliers[j]=getVar("multipliers",j);
    addMooseVariableDependency(_physical[j]);addMooseVariableDependency(_multipliers[j]);
    if(_physical[j]==_multipliers[j])paramError("multipliers","Algebraic variables must be distinct");
  }
  const auto allowed=getParam<std::vector<SubdomainID>>("allowed_primary_block_ids");
  if(allowed.empty())paramError("allowed_primary_block_ids","Source-approved host blocks required");
  std::ifstream input(getParam<FileName>("source_map"));
  if(!input)paramError("source_map","Cannot read source map");
  std::string line;std::getline(input,line);if(!line.empty() && line.back()=='\r')line.pop_back();
  if(line!="secondary_node,primary_element,secondary_block,x_m,y_m,z_m,host_node_0,host_node_1,host_node_2,host_node_3,host_node_4,host_node_5,host_node_6,host_node_7")
    paramError("source_map","Invalid schema");
  std::set<dof_id_type> seen;
  while(std::getline(input,line))
  {
    if(line.empty())continue;std::replace(line.begin(),line.end(),',',' ');std::istringstream values(line);
    Row row;dof_id_type eid,owner;Point expected;
    if(!(values>>row.nodes[0]>>eid>>owner>>expected(0)>>expected(1)>>expected(2)))paramError("source_map","Malformed row");
    for(unsigned j=1;j<9;++j)if(!(values>>row.nodes[j]))paramError("source_map","Missing host topology");
    std::string extra;if(values>>extra)paramError("source_map","Unexpected fields");
    if(!seen.insert(row.nodes[0]).second)paramError("source_map","Duplicate embedded node");
    auto node=_mesh.getMesh().query_node_ptr(row.nodes[0]);auto host=_mesh.getMesh().query_elem_ptr(eid);
    if(!node || !host || host->dim()!=3 || host->n_nodes()!=8)paramError("source_map","Replicated source node/HEX8 host missing");
    for(unsigned j=0;j<3;++j)
      if(!std::isfinite(expected(j)) || std::abs((*node)(j)-expected(j))>1e-11)paramError("source_map","Secondary coordinates differ");
    if(std::find(allowed.begin(),allowed.end(),host->subdomain_id())==allowed.end())paramError("source_map","Host outside declared blocks");
    for(unsigned j=0;j<8;++j)if(host->node_id(j)!=row.nodes[j+1])paramError("source_map","Host topology differs");
    auto it=_mesh.nodeToElemMap().find(row.nodes[0]);bool belongs=false;
    if(it!=_mesh.nodeToElemMap().end())for(auto child:it->second)
    {
      auto elem=_mesh.elemPtr(child);if(elem->subdomain_id()==owner && elem->dim()==1 && elem->n_nodes()==2)belongs=true;
    }
    if(!belongs || !host->contains_point(*node,1e-8))paramError("source_map","Invalid secondary owner or point outside source host");
    auto natural=libMesh::FEInterface::inverse_map(3,libMesh::FEType(FIRST,LAGRANGE),host,*node,1e-12,true);
    Point reconstructed;Real sum=0.;row.coefficients[0]=1.;
    for(unsigned j=0;j<8;++j)
    {
      Real weight=libMesh::FEInterface::shape(libMesh::FEType(FIRST,LAGRANGE),host,j,natural,false);
      row.coefficients[j+1]=-weight;sum+=weight;reconstructed+=weight*host->point(j);
    }
    if(std::abs(sum-1)>1e-10 || (reconstructed-*node).norm()>1e-8)paramError("source_map","Host shape interpolation identity failed");
    _subproblem.addGhostedElem(eid);
    if(node->processor_id()==processor_id())_rows.push_back(row);
  }
  if(seen.empty())paramError("source_map","Empty source map");
  // NodalConstraint invokes this batched assembler once. Avoid reinitializing
  // dense nodal arrays for every node in the entire reinforcement network.
  if(!_rows.empty())
  {
    _connected_nodes.push_back(_rows.front().nodes[0]);
    _primary_node_vector.push_back(_rows.front().nodes[1]);
  }
}
std::vector<dof_id_type> AqueductEmbeddedMPC::dofs(const Row & row,unsigned j) const
{
  std::vector<dof_id_type> ids;
  for(auto id:row.nodes)
  {
    auto node=_mesh.getMesh().query_node_ptr(id);
    if(!node->n_dofs(_sys.number(),_physical[j]->number()))mooseError("Embedded physical DOF missing");
    ids.push_back(node->dof_number(_sys.number(),_physical[j]->number(),0));
  }
  auto node=_mesh.getMesh().query_node_ptr(row.nodes[0]);
  if(!node->n_dofs(_sys.number(),_multipliers[j]->number()))mooseError("Embedded multiplier DOF missing");
  ids.push_back(node->dof_number(_sys.number(),_multipliers[j]->number(),0));return ids;
}
void AqueductEmbeddedMPC::computeResidual(const NumericVector<Number> &)
{
  const auto & solution=*_sys.currentSolution();
  for(const auto & row:_rows)for(unsigned j=0;j<3;++j)
  {
    const auto ids=dofs(row,j);DenseVector<Number> r(10);r.zero();Real gap=0;
    for(unsigned k=0;k<9;++k){r(k)=row.coefficients[k]*solution(ids[9])*_physical[j]->scalingFactor();gap+=row.coefficients[k]*solution(ids[k]);}
    r(9)=gap*_multipliers[j]->scalingFactor();addResiduals(_assembly,r,ids,1.);
  }
}
void AqueductEmbeddedMPC::computeJacobian(const SparseMatrix<Number> &)
{
  for(const auto & row:_rows)for(unsigned j=0;j<3;++j)
  {
    const auto ids=dofs(row,j);DenseMatrix<Number> K(10,10);K.zero();
    for(unsigned k=0;k<9;++k){K(k,9)=row.coefficients[k]*_physical[j]->scalingFactor();K(9,k)=row.coefficients[k]*_multipliers[j]->scalingFactor();}
    addJacobian(_assembly,K,ids,ids,1.);
  }
}
Real AqueductEmbeddedMPC::computeQpResidual(Moose::ConstraintType){mooseError("Batched embedded MPC uses full conservative rows");}
Real AqueductEmbeddedMPC::computeQpJacobian(Moose::ConstraintJacobianType){mooseError("Batched embedded MPC uses full conservative rows");}
