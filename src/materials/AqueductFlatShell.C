#include "AqueductFlatShell.h"
#include "libmesh/elem.h"

registerMooseObject("DamSafetyApp", AqueductFlatShellStress);
registerMooseObject("DamSafetyApp", AqueductFlatShell);
registerMooseObject("DamSafetyApp", AqueductFlatShellInertia);

namespace
{
AqueductShell::Quad quad(const Elem & elem)
{
  if(elem.n_nodes()!=4 || elem.dim()!=2)
    mooseError("AqueductFlatShell requires four-node rectangular shell elements");
  AqueductShell::Quad q;
  for(unsigned i=0;i<4;++i)
    for(unsigned j=0;j<3;++j) q[i][j]=elem.point(i)(j);
  return q;
}
AqueductShell::Vec position(const Point & p) {return {p(0),p(1),p(2)};}
void couplings(InputParameters & p)
{
  p.addRequiredCoupledVar("displacements","Three global displacement variables");
  p.addRequiredCoupledVar("rotations","Three global rotation variables, including drill");
  p.set<bool>("use_displaced_mesh")=false;
}
void section(InputParameters & p)
{
  p.addRequiredRangeCheckedParam<Real>("thickness","thickness>0","Physical shell thickness");
  p.addRequiredParam<Real>("offset","Reference-surface position / thickness; SPOS=0.5");
}
}

InputParameters AqueductFlatShellStress::validParams()
{
  auto p=Material::validParams();couplings(p);section(p);
  p.addClassDescription("Small-strain rectangular six-DOF MITC4 steel shell with offset and three Simpson layers");
  p.addRequiredRangeCheckedParam<Real>("youngs_modulus","youngs_modulus>0","Young modulus");
  p.addRequiredRangeCheckedParam<Real>("poissons_ratio","poissons_ratio>-1 & poissons_ratio<0.5","Poisson ratio");
  p.addRequiredRangeCheckedParam<Real>("yield_stress","yield_stress>0","Initial yield stress");
  p.addRequiredRangeCheckedParam<Real>("ultimate_stress","ultimate_stress>0","Last table stress");
  p.addRequiredRangeCheckedParam<Real>("plastic_strain_at_ultimate","plastic_strain_at_ultimate>0","Last table abscissa");
  p.addRequiredRangeCheckedParam<Real>("shear_correction","shear_correction>0 & shear_correction<=1","Explicit transverse shear correction");
  p.addRequiredRangeCheckedParam<Real>("drilling_factor","drilling_factor>0","Drilling stiffness / (G thickness); validate sensitivity separately");
  return p;
}
AqueductFlatShellStress::AqueductFlatShellStress(const InputParameters & p)
 :Material(p),_h(getParam<Real>("thickness")),_offset(getParam<Real>("offset")),
 _E(getParam<Real>("youngs_modulus")),_nu(getParam<Real>("poissons_ratio")),
 _fy(getParam<Real>("yield_stress")),_fu(getParam<Real>("ultimate_stress")),
 _limit(getParam<Real>("plastic_strain_at_ultimate")),_shear_factor(getParam<Real>("shear_correction")),
 _drill_factor(getParam<Real>("drilling_factor")),
 _state(declareProperty<std::vector<Real>>("aq_shell_state")),
 _old_state(getMaterialPropertyOld<std::vector<Real>>("aq_shell_state")),
 _strain(declareProperty<std::vector<Real>>("aq_shell_strain")),
 _old_strain(getMaterialPropertyOld<std::vector<Real>>("aq_shell_strain")),
 _membrane(declareADProperty<RealVectorValue>("aq_shell_membrane")),
 _moment(declareADProperty<RealVectorValue>("aq_shell_moment")),
 _shear(declareADProperty<RealVectorValue>("aq_shell_shear")),
 _drill(declareADProperty<Real>("aq_shell_drill"))
{
  if(coupledComponents("displacements")!=3 || coupledComponents("rotations")!=3)
    paramError("rotations","AqueductFlatShell needs three displacements and three rotations");
  if(_fu<_fy) paramError("ultimate_stress","Last table stress must not be less than initial yield");
  if(getParam<bool>("use_displaced_mesh")) paramError("use_displaced_mesh","SMALL strain requires the reference mesh");
  for(unsigned j=0;j<3;++j)
  {
    _dofs[j]=&adCoupledDofValues("displacements",j);
    _dofs[j+3]=&adCoupledDofValues("rotations",j);
  }
}
void AqueductFlatShellStress::initQpStatefulProperties()
{
  _state[_qp].assign(15,0.);_strain[_qp].assign(6,0.);
}
void AqueductFlatShellStress::computeQpProperties()
{
  auto B=AqueductShell::bMatrix(quad(*_current_elem),position(_q_point[_qp]));
  std::array<ADReal,9> g{};
  for(unsigned k=0;k<9;++k)
    for(unsigned i=0;i<4;++i)
      for(unsigned j=0;j<6;++j) g[k]+=B[k][6*i+j]*(*_dofs[j])[i];
  _state[_qp]=_old_state[_qp];_strain[_qp].resize(6);
  for(unsigned k=0;k<6;++k) _strain[_qp][k]=MetaPhysicL::raw_value(g[k]);
  _membrane[_qp].zero();_moment[_qp].zero();_shear[_qp].zero();
  auto z=AqueductShell::layerPositions(_h,_offset);
  for(unsigned t=0;t<3;++t)
  {
    std::array<ADReal,3> de;
    for(unsigned j=0;j<3;++j)
      de[j]=g[j]-_old_strain[_qp][j]+z[t]*(g[j+3]-_old_strain[_qp][j+3]);
    de[2]*=.5;
    std::array<double,4> old;
    for(unsigned j=0;j<4;++j) old[j]=_old_state[_qp][5*t+j];
    auto s=AqueductShell::planeStress(de,old,_old_state[_qp][5*t+4],_E,_nu,_fy,_fu,_limit);
    for(unsigned j=0;j<4;++j) _state[_qp][5*t+j]=MetaPhysicL::raw_value(s.value[j]);
    _state[_qp][5*t+4]=MetaPhysicL::raw_value(s.plastic_strain);
    for(unsigned j=0;j<3;++j)
    {
      const auto stress=s.value[j==2?3:j];
      _membrane[_qp](j)+=_h*AqueductShell::layerWeights[t]*stress;
      _moment[_qp](j)+=_h*AqueductShell::layerWeights[t]*z[t]*stress;
    }
  }
  const Real G=_E/(2*(1+_nu));
  for(unsigned j=0;j<2;++j) _shear[_qp](j)=_h*_shear_factor*G*g[j+6];
  _drill[_qp]=_h*_drill_factor*G*g[8];
}

InputParameters AqueductFlatShell::validParams()
{
  auto p=ADKernel::validParams();couplings(p);
  p.addRequiredRangeCheckedParam<unsigned int>("component","component<6","Global ux,uy,uz,rx,ry,rz component");
  return p;
}
AqueductFlatShell::AqueductFlatShell(const InputParameters & p)
 :ADKernel(p),_component(getParam<unsigned int>("component")),
 _membrane(getADMaterialProperty<RealVectorValue>("aq_shell_membrane")),
 _moment(getADMaterialProperty<RealVectorValue>("aq_shell_moment")),
 _shear(getADMaterialProperty<RealVectorValue>("aq_shell_shear")),
 _drill(getADMaterialProperty<Real>("aq_shell_drill"))
{
  if(getParam<bool>("use_displaced_mesh")) paramError("use_displaced_mesh","SMALL strain requires the reference mesh");
}
ADReal AqueductFlatShell::computeQpResidual()
{
  auto B=AqueductShell::bMatrix(quad(*_current_elem),position(_q_point[_qp]));
  unsigned col=6*_i+_component;
  ADReal residual=B[8][col]*_drill[_qp];
  for(unsigned j=0;j<3;++j)
    residual+=B[j][col]*_membrane[_qp](j)+B[j+3][col]*_moment[_qp](j);
  for(unsigned j=0;j<2;++j) residual+=B[j+6][col]*_shear[_qp](j);
  return residual;
}

InputParameters AqueductFlatShellInertia::validParams()
{
  auto p=ADTimeKernel::validParams();couplings(p);section(p);
  p.addRequiredRangeCheckedParam<Real>("density","density>0","Mass density");
  p.addRequiredRangeCheckedParam<unsigned int>("component","component<6","Global displacement/rotation component");
  return p;
}
AqueductFlatShellInertia::AqueductFlatShellInertia(const InputParameters & p)
 :ADTimeKernel(p),_component(getParam<unsigned int>("component")),_rho(getParam<Real>("density")),
 _mass(AqueductShell::moments(getParam<Real>("thickness"),getParam<Real>("offset")))
{
  if(coupledComponents("displacements")!=3 || coupledComponents("rotations")!=3)
    paramError("rotations","AqueductFlatShell needs three displacements and three rotations");
  if(getParam<bool>("use_displaced_mesh")) paramError("use_displaced_mesh","SMALL strain requires the reference mesh");
  for(unsigned j=0;j<3;++j)
  {
    _accelerations[j]=&adCoupledDotDot("displacements",j);
    _accelerations[j+3]=&adCoupledDotDot("rotations",j);
  }
}
ADReal AqueductFlatShellInertia::computeQpResidual()
{
  auto n=AqueductShell::frame(quad(*_current_elem)).n;
  std::array<ADReal,3> a,r;
  ADReal normal_rotation=0;
  for(unsigned j=0;j<3;++j)
  {
    a[j]=(*_accelerations[j])[_qp];r[j]=(*_accelerations[j+3])[_qp];
    normal_rotation+=n[j]*r[j];
  }
  unsigned j=_component%3,k=(j+1)%3,l=(j+2)%3;
  const ADReal residual=_component<3 ?
    _mass.zero*a[j]+_mass.first*(r[k]*n[l]-r[l]*n[k]) :
    _mass.first*(n[k]*a[l]-n[l]*a[k])+_mass.second*(r[j]-n[j]*normal_rotation);
  return _rho*_test[_i][_qp]*residual;
}
