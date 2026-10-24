#include "AqueductHHT.h"
#include "AqueductHHTScheme.h"
#include "NonlinearSystem.h"
#include "FEProblemBase.h"
#include "libmesh/dof_map.h"
#include "libmesh/nonlinear_implicit_system.h"

registerMooseObject("DamSafetyApp", AqueductHHT);
InputParameters AqueductHHT::validParams()
{
  auto p=NewmarkBeta::validParams();
  p.addRequiredRangeCheckedParam<Real>("hht_alpha","hht_alpha>=-0.5 & hht_alpha<=0","Explicit HHT alpha; no assumed Abaqus version default");
  p.addParam<Real>("hht_start_time",0.,"Unweighted quasi-static stage ends at this total time");
  p.addParam<std::vector<VariableName>>("algebraic_variables",{},"Multiplier equations stay at the current time, with zero old-force entries");
  p.addClassDescription("Whole-system normalized HHT accepted-force history; pair every mass with HHT-scaled inertia");
  return p;
}
AqueductHHT::AqueductHHT(const InputParameters & p)
 :NewmarkBeta(p),_hht_alpha(getParam<Real>("hht_alpha")),_start_time(getParam<Real>("hht_start_time")),
 _old_forces(addVector("aq_hht_accepted_forces",true,libMesh::GHOSTED))
{
  if(_var_restriction) paramError("variables","Whole-system HHT cannot restrict physical variables");
  if(std::abs(_beta-AqueductHHTScheme::beta(_hht_alpha))>1e-12 ||
     std::abs(_gamma-AqueductHHTScheme::gamma(_hht_alpha))>1e-12)
    paramError("beta","Set beta=(1-alpha)^2/4 and gamma=.5-alpha consistently in integrator, BCs and auxiliaries");
}
void AqueductHHT::initialSetup()
{
  _algebraic_dofs.clear();
  if(!_nl)return; // MOOSE also constructs an auxiliary-system integrator.
  auto & system=_sys.system();
  for(const auto & name:getParam<std::vector<VariableName>>("algebraic_variables"))
  {
    if(!system.has_variable(name))paramError("algebraic_variables","Unknown multiplier variable "+name);
    std::vector<dof_id_type> ids;
    system.get_dof_map().local_variable_indices(ids,system.get_mesh(),system.variable_number(name));
    _algebraic_dofs.insert(_algebraic_dofs.end(),ids.begin(),ids.end());
  }
}
void AqueductHHT::clearAlgebraicHistory()
{
  for(auto dof:_algebraic_dofs)_old_forces->set(dof,0.);
  _old_forces->close();
}
void AqueductHHT::init()
{
  NewmarkBeta::init();_old_forces->zero();_old_forces->close();
  if(_nl)
  {
    // Same initial non-time assembly contract as locked CrankNicolson.
    _fe_problem.setCurrentResidualVectorTags({_nl->nonTimeVectorTag()});
    _nl->computeResidualTag(_nl->RHS(),_nl->nonTimeVectorTag());
    _fe_problem.clearCurrentResidualVectorTags();copyVector(_nl->RHS(),*_old_forces);
    clearAlgebraicHistory();
  }
}
void AqueductHHT::postResidual(NumericVector<Number> & residual)
{
  std::vector<unsigned char> closed={_Re_time->closed(),_Re_non_time->closed(),_old_forces->closed()};
  comm().min(closed);
  if(!closed[0])_Re_time->close();
  if(!closed[1])_Re_non_time->close();
  if(!closed[2])_old_forces->close();
  residual+=*_Re_time;residual+=*_Re_non_time;
  if(_fe_problem.time()>_start_time)
    residual.add(-_hht_alpha*AqueductHHTScheme::massScale(_hht_alpha),*_old_forces);
  residual.close();
}
void AqueductHHT::postStep()
{
  // TransientBase calls this only for an accepted solve, never on a rejected attempt.
  // Contact-search EquationSystems::reinit must preserve this accepted physical history.
  copyVector(*_Re_non_time,*_old_forces);clearAlgebraicHistory();
}
void AqueductHHT::computeADTimeDerivatives(ADReal & v,const dof_id_type & dof,ADReal & a) const
{
  if(_t_step<=_inactive_tsteps){v=0.;a=0.;}
  else NewmarkBeta::computeADTimeDerivatives(v,dof,a);
}
