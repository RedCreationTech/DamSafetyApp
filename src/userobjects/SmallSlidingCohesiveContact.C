#include "SmallSlidingCohesiveContact.h"
#include "MooseVariable.h"
#include "SystemBase.h"
#include "Assembly.h"
#include "AutomaticMortarGeneration.h"
#include "MortarContactUtils.h"
#include "libmesh/quadrature.h"

registerMooseObject("DamSafetyApp", SmallSlidingCohesiveContact);

InputParameters SmallSlidingCohesiveContact::validParams()
{
  auto p=MortarUserObject::validParams();
  p.addClassDescription("Reference-geometry small-sliding MAXS/displacement cohesive contact.");
  p.addRequiredCoupledVar("displacements","Three displacement fields.");
  p.addRequiredCoupledVar("normal_pressure","Normal HARD-contact multiplier in Pa.");
  p.addRequiredParam<std::vector<Real>>("stiffness","Kn Ks Kt, in Pa/m.");
  p.addRequiredParam<std::vector<Real>>("strength","Three MAXS strengths, in Pa.");
  p.addRequiredRangeCheckedParam<Real>("failure_increment","failure_increment > 0","Final minus onset effective separation, m.");
  p.addRangeCheckedParam<Real>("viscosity",0,"viscosity >= 0","Damage stabilization time, s.");
  p.addRangeCheckedParam<Real>("friction_coefficient",0,"friction_coefficient >= 0","Original Coulomb coefficient.");
  p.addRequiredRangeCheckedParam<Real>("critical_elastic_slip","critical_elastic_slip > 0","Slip tolerance times reference contact-element length, m.");
  p.addRangeCheckedParam<Real>("initial_bond_tolerance",1e-10,"initial_bond_tolerance > 0","Geometry roundoff tolerance for original-contact eligibility, m.");
  p.addParam<unsigned int>("expected_secondary_nodes",0,"Strict full interface coverage count; zero disables the count check.");
  p.set<bool>("use_displaced_mesh")=false;
  p.set<bool>("interpolate_normals")=false;
  p.set<ExecFlagEnum>("execute_on")={EXEC_LINEAR,EXEC_NONLINEAR,EXEC_TIMESTEP_END};
  return p;
}

SmallSlidingCohesiveContact::SmallSlidingCohesiveContact(const InputParameters & p)
 : MortarUserObject(p),
   _secondary{&adCoupledValue("displacements",0),&adCoupledValue("displacements",1),&adCoupledValue("displacements",2)},
   _primary{&adCoupledNeighborValue("displacements",0),&adCoupledNeighborValue("displacements",1),&adCoupledNeighborValue("displacements",2)},
   _displacement_vars{getVar("displacements",0),getVar("displacements",1),getVar("displacements",2)},
   _pressure_var(getVar("normal_pressure",0)),
   _system(*getCheckedPointerParam<SystemBase *>("_sys")),
   _coord(_mci_assembly.mortarCoordTransformation()),
   _parameters{{getParam<std::vector<Real>>("stiffness").at(0),getParam<std::vector<Real>>("stiffness").at(1),getParam<std::vector<Real>>("stiffness").at(2)},
               {getParam<std::vector<Real>>("strength").at(0),getParam<std::vector<Real>>("strength").at(1),getParam<std::vector<Real>>("strength").at(2)},
               getParam<Real>("failure_increment"),getParam<Real>("viscosity"),getParam<Real>("friction_coefficient"),getParam<Real>("critical_elastic_slip")},
   _bond_tolerance(getParam<Real>("initial_bond_tolerance")),
   _expected_nodes(getParam<unsigned int>("expected_secondary_nodes")),
   _history(declareRestartableData<std::unordered_map<dof_id_type,AbaqusCohesiveLaw::State>>("accepted_contact_history"))
{
  if (getParam<bool>("use_displaced_mesh") || interpolateNormals())
    mooseError("SmallSlidingCohesiveContact requires fixed reference projection/nodal normals.");
  if(coupledComponents("displacements")!=3 || getParam<std::vector<Real>>("stiffness").size()!=3 || getParam<std::vector<Real>>("strength").size()!=3)
    mooseError("Contact requires three displacement, stiffness and strength components.");
  for (unsigned int d=0;d<3;++d)
    if(_parameters.stiffness[d]<=0 || _parameters.strength[d]<=0)
      mooseError("Cohesive stiffnesses and MAXS strengths must be positive.");
}
void SmallSlidingCohesiveContact::initialSetup()
{
  MortarUserObject::initialSetup();
  _shape=&_displacement_vars[0]->phiLower();
}
void SmallSlidingCohesiveContact::initialize()
{
  _gaps.clear(); _jumps.clear(); _reference_normals.clear(); _initial_gaps.clear();
}
void SmallSlidingCohesiveContact::execute()
{
  const auto & primary_map=amg().getPrimaryIpToLowerElementMap(*_lower_primary_elem,*_lower_primary_elem->interior_parent(),*_lower_secondary_elem);
  const auto & secondary_map=amg().getSecondaryIpToLowerElementMap(*_lower_secondary_elem);
  for(unsigned int qp=0;qp<_qrule_msm->n_points();++qp)
  {
    std::array<ADReal,3> prim{(*_primary[0])[qp],(*_primary[1])[qp],(*_primary[2])[qp]};
    std::array<ADReal,3> sec{(*_secondary[0])[qp],(*_secondary[1])[qp],(*_secondary[2])[qp]};
    trimInteriorNodeDerivatives(primary_map,_displacement_vars,prim,false);
    trimInteriorNodeDerivatives(secondary_map,_displacement_vars,sec,true);
    const ADRealVectorValue jump(prim[0]-sec[0],prim[1]-sec[1],prim[2]-sec[2]);
    const Point initial=_phys_points_primary[qp]-_phys_points_secondary[qp];
    for(unsigned int i=0;i<_shape->size();++i)
    {
      const auto * node=_lower_secondary_elem->node_ptr(i);
      const Real weight=(*_shape)[i][qp]*_JxW_msm[qp]*_coord[qp];
      _gaps[node].first+=weight*(jump*_normals[i]+initial*_normals[i]);
      _gaps[node].second+=weight;
      _initial_gaps[node]+=weight*(initial*_normals[i]);
      _jumps[node]+=weight*jump;
      _reference_normals[node]+=weight*ADRealVectorValue(_normals[i]);
    }
  }
}
void SmallSlidingCohesiveContact::finalize()
{
  using namespace Moose::Mortar::Contact;
  communicateGaps(_gaps,_mci_mesh,true,true,_communicator,true);
  communicateRealObject(_jumps,_mci_mesh,true,_communicator,true);
  communicateRealObject(_reference_normals,_mci_mesh,true,_communicator,true);
  communicateRealObject(_initial_gaps,_mci_mesh,true,_communicator,true);
  if(!_coverage_reported)
  {
    unsigned int covered=0,bonded=0;Real area=0;
    for(const auto & entry:_gaps)
      if(entry.first->processor_id()==processor_id())
      {
        ++covered;area+=entry.second.second;
        if(_initial_gaps.at(entry.first)/entry.second.second<=_bond_tolerance)++bonded;
      }
    _communicator.sum(covered);_communicator.sum(bonded);_communicator.sum(area);
    if(_expected_nodes && covered!=_expected_nodes)
      mooseError(name()," reference mortar coverage is ",covered," nodes, expected ",_expected_nodes);
    mooseInfo(name()," initial interface coverage: ",covered," nodes, original-contact eligible: ",bonded,"; area m2: ",area);
    _coverage_reported=true;
  }
  if(_mci_fe_problem.getCurrentExecuteOnFlag()==EXEC_TIMESTEP_END)
    for(const auto & entry : _gaps)
      _history[entry.first->id()]=evaluateNode(static_cast<const Node *>(entry.first)).state;
}
std::array<Point,3> SmallSlidingCohesiveContact::basis(const Node * node) const
{
  Point normal=MetaPhysicL::raw_value(_reference_normals.at(node)); normal/=normal.norm();
  Point tangent(1,0,0); tangent-=normal(0)*normal;
  if(tangent.norm()<1e-8) { tangent=Point(0,1,0);tangent-=normal(1)*normal; }
  tangent/=tangent.norm();
  return {normal,tangent,normal.cross(tangent)};
}
AbaqusCohesiveLaw::Result<ADReal> SmallSlidingCohesiveContact::evaluateNode(const Node * node)
{
  const auto rotation=basis(node);
  const auto jump=_jumps.at(node)/_gaps.at(node).second;
  std::array<ADReal,3> local{jump*rotation[0],jump*rotation[1],jump*rotation[2]};
  const auto dof=node->dof_number(_system.number(),_pressure_var->number(),0);
  ADReal pressure=(*_system.currentSolution())(dof); Moose::derivInsert(pressure.derivatives(),dof,1.0);
  const bool bonded=_initial_gaps.at(node)/_gaps.at(node).second<=_bond_tolerance;
  return AbaqusCohesiveLaw::evaluate(_parameters,_history[node->id()],local,pressure,_mci_fe_problem.dt(),bonded);
}
void SmallSlidingCohesiveContact::reinit()
{
  for(auto & field:_traction)
  {
    field.resize(_qrule_msm->n_points());
    for(unsigned int qp=0;qp<field.size();++qp)field[qp]=0;
  }
  for(unsigned int i=0;i<_shape->size();++i)
  {
    const auto * node=_lower_secondary_elem->node_ptr(i);
    if(!_gaps.count(node)) continue;
    const auto local=evaluateNode(node); const auto frame=basis(node);
    for(unsigned int qp=0;qp<_qrule_msm->n_points();++qp)
      for(unsigned int d=0;d<3;++d)
        for(unsigned int k=0;k<3;++k)
          _traction[d][qp]+=(*_shape)[i][qp]*frame[k](d)*local.traction[k];
  }
}
Real SmallSlidingCohesiveContact::nodeValue(const Node * node,unsigned int quantity) const
{
  const auto it=_gaps.find(node);
  if(it==_gaps.end())return 0;
  const auto history=_history.find(node->id());
  const Real damage=history==_history.end() ? 0:history->second[7];
  if(quantity==0)return damage;
  if(quantity==1)return MetaPhysicL::raw_value(it->second.first)/it->second.second;
  const auto dof=node->dof_number(_system.number(),_pressure_var->number(),0);
  const Real pressure=(*_system.currentSolution())(dof);
  if(quantity==2)return pressure;
  const auto frame=basis(node);const auto jump=MetaPhysicL::raw_value(_jumps.at(node))/it->second.second;
  if(quantity>=3 && quantity<=5)return jump*frame[quantity-3];
  const unsigned int d=quantity-6;
  if(d>2)mooseError("Unknown contact output quantity.");
  const Real separation=jump*frame[d];
  Real residual=-(1-damage)*_parameters.stiffness[d]*(d==0 ? std::max(separation,0.0):separation);
  if(d==0)residual+=pressure;
  else if(history!=_history.end())residual-=damage*_parameters.friction*std::max(pressure,0.0)*history->second[4+d]/_parameters.critical_slip;
  return residual;
}
