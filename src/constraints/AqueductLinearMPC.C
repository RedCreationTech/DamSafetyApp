#include "AqueductLinearMPC.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SystemBase.h"
#include "Assembly.h"
#include "libmesh/elem.h"
#include "libmesh/sparse_matrix.h"
#include <set>
#include <cmath>

registerMooseObject("DamSafetyApp", AqueductLinearMPC);
registerMooseObject("DamSafetyApp", AqueductConstraintNodesGenerator);

InputParameters AqueductLinearMPC::validParams()
{
  auto p=NodalConstraint::validParams();
  p.addClassDescription("Conservative linear multivariable mixed row with optional physical compliance");
  p.addRequiredParam<dof_id_type>("anchor_node","Node carrying this row's multiplier");
  p.addRequiredParam<std::vector<dof_id_type>>("term_nodes","Node for every linear term");
  p.addRequiredCoupledVar("terms","Variable for every term, in the same order as term_nodes");
  p.addRequiredParam<std::vector<Real>>("coefficients","Coefficients in C*u=0");
  p.addRequiredCoupledVar("multiplier","One multiplier variable active at anchor_node");
  p.addRangeCheckedParam<Real>("multiplier_compliance",0.,"multiplier_compliance>=0","Physical compliance in C*u-compliance*lambda=0; zero is a kinematic row");
  p.addParam<bool>("reuse_assembly_storage",false,"Reuse constant row matrices and containers; refresh DOF values every call");
  p.suppressParameter<MooseEnum>("formulation");
  p.suppressParameter<NonlinearVariableName>("variable_secondary");
  p.set<bool>("use_displaced_mesh")=false;
  return p;
}
AqueductLinearMPC::AqueductLinearMPC(const InputParameters & p)
 :NodalConstraint(p),_anchor(getParam<dof_id_type>("anchor_node")),
 _nodes(getParam<std::vector<dof_id_type>>("term_nodes")),
 _coefficients(getParam<std::vector<Real>>("coefficients")),_multiplier(getVar("multiplier",0)),_compliance(getParam<Real>("multiplier_compliance"))
{
  if(_nodes.size()<2 || _nodes.size()!=_coefficients.size() || _nodes.size()!=coupledComponents("terms"))
    paramError("term_nodes","At least two matching nodes, variables and coefficients are required");
  if(coupledComponents("multiplier")!=1) paramError("multiplier","Exactly one multiplier is required");
  if(getParam<bool>("use_displaced_mesh")) paramError("use_displaced_mesh","Linear MPC uses reference coordinates");
  if(!std::isfinite(_compliance))paramError("multiplier_compliance","Finite physical compliance required");
  std::set<dof_id_type> primary;
  for(unsigned j=0;j<_nodes.size();++j)
  {
    if(!std::isfinite(_coefficients[j]))paramError("coefficients","Nonfinite coefficient");
    _variables.push_back(getVar("terms",j));addMooseVariableDependency(_variables.back());
    if(_variables.back()==_multiplier)paramError("terms","Multiplier cannot be a physical term");
    primary.insert(_nodes[j]);
  }
  if(primary.empty() || (primary.size()==1 && *primary.begin()==_anchor))
    paramError("term_nodes","At least one non-anchor node is required");
  auto & mesh=_mesh.getMesh();auto node=mesh.query_node_ptr(_anchor);
  if(!node)paramError("anchor_node","Anchor not present; replicated mesh is required");
  if(node->processor_id()==processor_id())_connected_nodes.push_back(_anchor);
  _primary_node_vector.assign(primary.begin(),primary.end());
  for(auto id:primary)
  {
    if(!mesh.query_node_ptr(id))paramError("term_nodes","Missing primary node");
    auto it=_mesh.nodeToElemMap().find(id);
    if(it!=_mesh.nodeToElemMap().end())for(auto eid:it->second)_subproblem.addGhostedElem(eid);
  }
  addMooseVariableDependency(_multiplier);
  if(getParam<bool>("reuse_assembly_storage"))
  {
    _storage=std::make_unique<AqueductMPCAssemblyStorage>(_coefficients);
    _scales_scratch.resize(_coefficients.size());
  }
}
std::vector<dof_id_type> AqueductLinearMPC::dofs()
{
  std::vector<dof_id_type> result;
  dofs(result);return result;
}
void AqueductLinearMPC::dofs(std::vector<dof_id_type> & result)
{
  result.clear();
  for(unsigned j=0;j<_nodes.size();++j)
  {
    auto node=_mesh.getMesh().query_node_ptr(_nodes[j]);
    if(!node || !node->n_dofs(_sys.number(),_variables[j]->number()))
      mooseError("AqueductLinearMPC term variable has no DOF on specified node");
    result.push_back(node->dof_number(_sys.number(),_variables[j]->number(),0));
  }
  auto node=_mesh.getMesh().query_node_ptr(_anchor);
  if(!node->n_dofs(_sys.number(),_multiplier->number()))
    mooseError("AqueductLinearMPC multiplier has no DOF on anchor");
  result.push_back(node->dof_number(_sys.number(),_multiplier->number(),0));
}
void AqueductLinearMPC::computeResidual(const NumericVector<Number> &)
{
  if(_storage)
  {
    auto & ids=_storage->dofs();dofs(ids);const auto & solution=*_sys.currentSolution();
    auto & r=_storage->residual();Real gap=0;
    for(unsigned j=0;j<_nodes.size();++j)
    {
      r(j)=_coefficients[j]*solution(ids.back())*_variables[j]->scalingFactor();
      gap+=_coefficients[j]*solution(ids[j]);
    }
    r(ids.size()-1)=(gap-_compliance*solution(ids.back()))*_multiplier->scalingFactor();
    addResiduals(_assembly,r,ids,1.);return;
  }
  auto ids=dofs();const auto & solution=*_sys.currentSolution();
  DenseVector<Number> r(ids.size());r.zero();Real gap=0;
  for(unsigned j=0;j<_nodes.size();++j)
  {
    r(j)=_coefficients[j]*solution(ids.back())*_variables[j]->scalingFactor();
    gap+=_coefficients[j]*solution(ids[j]);
  }
  r(ids.size()-1)=(gap-_compliance*solution(ids.back()))*_multiplier->scalingFactor();
  addResiduals(_assembly,r,ids,1.);
}
void AqueductLinearMPC::computeJacobian(const SparseMatrix<Number> &)
{
  if(_storage)
  {
    auto & ids=_storage->dofs();dofs(ids);
    for(unsigned j=0;j<_nodes.size();++j)_scales_scratch[j]=_variables[j]->scalingFactor();
    const auto & K=_storage->jacobian(_scales_scratch,_multiplier->scalingFactor(),_compliance);
    addJacobian(_assembly,K,ids,ids,1.);return;
  }
  auto ids=dofs();DenseMatrix<Number> K(ids.size(),ids.size());K.zero();
  for(unsigned j=0;j<_nodes.size();++j)
  {
    K(j,ids.size()-1)=_coefficients[j]*_variables[j]->scalingFactor();
    K(ids.size()-1,j)=_coefficients[j]*_multiplier->scalingFactor();
  }
  K(ids.size()-1,ids.size()-1)=-_compliance*_multiplier->scalingFactor();
  addJacobian(_assembly,K,ids,ids,1.);
}
Real AqueductLinearMPC::computeQpResidual(Moose::ConstraintType)
{mooseError("AqueductLinearMPC assembles the complete conservative nodal row");}
Real AqueductLinearMPC::computeQpJacobian(Moose::ConstraintJacobianType)
{mooseError("AqueductLinearMPC assembles the complete conservative nodal row");}

InputParameters AqueductConstraintNodesGenerator::validParams()
{
  auto p=MeshGenerator::validParams();
  p.addRequiredParam<MeshGeneratorName>("input","Input reference mesh");
  p.addRequiredParam<std::vector<dof_id_type>>("nodes","Existing anchor node IDs");
  p.addRequiredParam<SubdomainID>("block_id","Unused block ID for multiplier support");
  p.addRequiredParam<SubdomainName>("block_name","Multiplier-support block name");
  return p;
}
AqueductConstraintNodesGenerator::AqueductConstraintNodesGenerator(const InputParameters & p)
 :MeshGenerator(p),_input(getMesh("input")){}
std::unique_ptr<MeshBase> AqueductConstraintNodesGenerator::generate()
{
  auto mesh=std::move(_input);auto block=getParam<SubdomainID>("block_id");
  for(const auto * elem:mesh->active_element_ptr_range())
    if(elem->subdomain_id()==block)paramError("block_id","Block ID already used");
  std::set<dof_id_type> unique;
  for(auto id:getParam<std::vector<dof_id_type>>("nodes"))
  {
    if(!unique.insert(id).second)paramError("nodes","Duplicate anchor");
    auto node=mesh->query_node_ptr(id);if(!node)paramError("nodes","Missing anchor node");
    auto elem=Elem::build(libMesh::NODEELEM);elem->set_node(0)=node;
    elem->subdomain_id()=block;elem->processor_id()=node->processor_id();
    elem->set_id(mesh->max_elem_id());mesh->add_elem(std::move(elem));
  }
  mesh->subdomain_name(block)=getParam<SubdomainName>("block_name");
  return mesh;
}
