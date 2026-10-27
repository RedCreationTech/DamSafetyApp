#include "ADDS04GSP1Correction.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "MooseException.h"
#include "RankFourTensor.h"
#include "libmesh/elem.h"
#include <algorithm>

registerMooseObject("DamSafetyApp", ADDS04GSP1Correction);
InputParameters ADDS04GSP1Correction::validParams()
{
  auto p=ADKernel::validParams();
  p.addClassDescription("Default-off DS04G center rule correction, reference HG/mass and current gravity; no constitutive update.");
  p.addRequiredParam<MooseEnum>("term",MooseEnum("hourglass mass gravity"),"Compensated operator");
  p.addRequiredParam<unsigned>("component","Displacement component");
  p.addRequiredCoupledVar("displacements","Two absolute displacement variables");
  p.addRequiredCoupledVar("velocity","Newmark velocity");
  p.addRequiredCoupledVar("acceleration","Newmark acceleration");
  for(const auto & n:{"youngs_modulus","poissons_ratio","expected_density","eta","zeta","alpha","beta","gamma","static_end","gravity"})
    p.addRequiredParam<Real>(n,"Frozen DS04G constant");
  p.addParam<MaterialPropertyName>("density","density","Constant density");
  p.addParam<MaterialPropertyName>("elasticity_tensor","elasticity_tensor","Constant undamaged C0");
  p.set<bool>("use_displaced_mesh")=false;
  return p;
}
ADDS04GSP1Correction::ADDS04GSP1Correction(const InputParameters & p)
 : ADKernel(p),_term(getParam<MooseEnum>("term")),_component(getParam<unsigned>("component")),
 _youngs(getParam<Real>("youngs_modulus")),_nu(getParam<Real>("poissons_ratio")),
 _rho(getParam<Real>("expected_density")),_eta(getParam<Real>("eta")),_zeta(getParam<Real>("zeta")),
 _alpha(getParam<Real>("alpha")),_beta(getParam<Real>("beta")),_gamma(getParam<Real>("gamma")),
 _static_end(getParam<Real>("static_end")),_gravity(getParam<Real>("gravity")),
 _density(getMaterialProperty<Real>("density")),_elasticity(getMaterialProperty<RankFourTensor>("elasticity_tensor")),
 _old(_var.dofValuesOld()),_older(_var.dofValuesOlder()),
 _v_old(coupledDofValuesOld("velocity")),_a_old(coupledDofValuesOld("acceleration"))
{
  if (_mesh.dimension()!=2 || _component>1 || _use_displaced_mesh || coupledComponents("displacements")!=2 ||
      !(_youngs>0 && _rho>0 && _nu>-1 && _nu<.5 && _beta>0 && _gamma>0 && _static_end>=0 && _zeta>=0 && _eta>=0 && _alpha>=-1./3 && _alpha<=0) || !std::isfinite(_gravity))
    mooseError("SP1 requires 2D XY, reference assembly and valid constant parameters");
  if(_var.number()!=coupled("displacements",_component)) paramError("variable","Must match component");
  for(unsigned c=0;c<2;++c) _u[c]=&adCoupledDofValues("displacements",c);
}
void ADDS04GSP1Correction::initialSetup()
{
  if(getBlockCoordSystem()!=Moose::COORD_XYZ) mooseError("SP1 requires Cartesian XY");
  for(unsigned c=0;c<2;++c)
  {
    auto * u=getVar("displacements",c);
    if(u->feType().family!=libMesh::LAGRANGE || u->feType().order!=libMesh::FIRST)
      mooseError("SP1 requires FIRST LAGRANGE displacements");
  }
}
ADReal ADDS04GSP1Correction::computeQpResidual()
{
  if(_qrule->n_points()!=1 || _current_elem->type()!=QUAD4 || _u[0]->size()!=4 || _u[1]->size()!=4)
    mooseError("SP1 requires one center QP and exactly four-node QUAD4");
  if(!std::isfinite(_density[_qp]) || std::abs(_density[_qp]-_rho)>1e-12*_rho)
    mooseError("SP1 density differs from frozen constant");
  const Real mu=_youngs/(2*(1+_nu)), lam=_youngs*_nu/((1+_nu)*(1-2*_nu));
  for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j) for(unsigned k=0;k<3;++k) for(unsigned l=0;l<3;++l)
  {
    const Real v=lam*(i==j)*(k==l)+mu*((i==k)*(j==l)+(i==l)*(j==k));
    if(!std::isfinite(_elasticity[_qp](i,j,k,l)) || std::abs(_elasticity[_qp](i,j,k,l)-v)>1e-10*_youngs)
      mooseError("SP1 requires declared undamaged isotropic C0");
  }
  DS04GSP1::Nodes nodes;
  for(unsigned i=0;i<4;++i) for(unsigned c=0;c<2;++c) nodes[i][c]=_current_elem->point(i)(c);
  try
  {
    const auto g=DS04GSP1::geometry(nodes);
    if(_term=="gravity")
    {
      if(_component!=1) mooseError("SP1 gravity is Y only");
      std::array<std::array<ADReal,2>,4> current;
      for(unsigned i=0;i<4;++i) for(unsigned c=0;c<2;++c) current[i][c]=nodes[i][c]+(*_u[c])[i];
      return -_rho*_gravity*DS04GSP1::gravityCorrection(current)[_i]/g.area;
    }
    const bool dynamic=_t>_static_end;
    if(_term=="mass")
    {
      if(!dynamic || _dt<=0) return 0.;
      ADReal f=0.;
      for(unsigned j=0;j<4;++j)
      {
        const ADReal a=DS04GSP1::trialAcceleration((*_u[_component])[j],_old[j],_v_old[j],_a_old[j],_dt,_beta);
        const ADReal v=_v_old[j]+_dt*(1-_gamma)*_a_old[j]+_gamma*_dt*a;
        f+=g.mass_correction[_i][j]*(a+_eta*((1+_alpha)*v-_alpha*_v_old[j]));
      }
      return _rho*f/g.area;
    }
    std::array<ADReal,4> u;std::array<Real,4> old,older;
    for(unsigned j=0;j<4;++j) {u[j]=(*_u[_component])[j];old[j]=_old[j];older[j]=_older[j];}
    ADReal f=DS04GSP1::hourglass(g,u,_i,mu);
    if(dynamic)
    {
      const auto w=DS04GSP1::history(_alpha,_zeta,_dt);
      f=w.current*f+w.old*DS04GSP1::hourglass(g,old,_i,mu)+w.older*DS04GSP1::hourglass(g,older,_i,mu);
    }
    return f/g.area;
  }
  catch(const std::invalid_argument & e) {throw MooseException(e.what());}
}
