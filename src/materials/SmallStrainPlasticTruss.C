#include "SmallStrainPlasticTruss.h"
#include "MooseVariable.h"
#include "libmesh/quadrature.h"
registerMooseObject("DamSafetyApp", SmallStrainPlasticTruss);
InputParameters SmallStrainPlasticTruss::validParams()
{
  auto p=TrussMaterial::validParams();
  p.addClassDescription("Reference T3D2 axial strain and two-point isotropic hardening with final plateau.");
  p.addRequiredRangeCheckedParam<Real>("yield_stress","yield_stress > 0","Initial yield stress, Pa.");
  p.addRequiredRangeCheckedParam<Real>("ultimate_stress","ultimate_stress > 0","Final tabulated stress, Pa.");
  p.addRequiredRangeCheckedParam<Real>("plastic_strain_knot","plastic_strain_knot > 0","Second tabulated equivalent plastic strain.");
  p.set<bool>("use_displaced_mesh")=false;
  return p;
}
SmallStrainPlasticTruss::SmallStrainPlasticTruss(const InputParameters & p)
 : TrussMaterial(p),_yield(getParam<Real>("yield_stress")),_ultimate(getParam<Real>("ultimate_stress")),_knot(getParam<Real>("plastic_strain_knot")),
   _old_strain(getMaterialPropertyOld<Real>(_base_name+"total_stretch")),
   _old_stress(getMaterialPropertyOld<Real>(_base_name+"axial_stress")),
   _plastic(declareProperty<Real>(_base_name+"plastic_stretch")),
   _equivalent(declareProperty<Real>(_base_name+"equivalent_plastic_strain")),
   _old_plastic(getMaterialPropertyOld<Real>(_base_name+"plastic_stretch")),
   _old_equivalent(getMaterialPropertyOld<Real>(_base_name+"equivalent_plastic_strain"))
{
  if (_ndisp!=3 || _ultimate<_yield || getParam<bool>("use_displaced_mesh"))
    mooseError("Small strain truss requires three reference displacements and nondecreasing stress table.");
}
void SmallStrainPlasticTruss::initQpStatefulProperties()
{
  TrussMaterial::initQpStatefulProperties();_plastic[_qp]=0;_equivalent[_qp]=0;
}
void SmallStrainPlasticTruss::computeProperties()
{
  if (_current_elem->n_nodes()!=2) mooseError("T3D2 requires two nodes.");
  const Point direction=_current_elem->point(1)-_current_elem->point(0);
  const Real length=direction.norm();Real strain=0;
  for(unsigned int d=0;d<3;++d)
    strain+=direction(d)*(_disp_var[d]->getNodalValue(_current_elem->node_ref(1))-
                         _disp_var[d]->getNodalValue(_current_elem->node_ref(0)))/(length*length);
  const Real H=(_ultimate-_yield)/_knot;
  for(_qp=0;_qp<_qrule->n_points();++_qp)
  {
    const Real E=_youngs_modulus[_qp];
    const Real trial=_old_stress[_qp]+E*(strain-_old_strain[_qp]);
    const Real old_p=_old_equivalent[_qp];
    const Real f=std::abs(trial)-(_yield+H*std::min(old_p,_knot));
    Real dp=0,tangent=E;
    if(f>0)
    {
      dp=f/(E+(old_p<_knot ? H:0));
      if(old_p+dp>=_knot) {dp=(std::abs(trial)-_ultimate)/E;tangent=0;}
      else tangent=E*H/(E+H);
    }
    const Real sign=trial>=0 ? 1:-1;
    _total_stretch[_qp]=strain;
    _equivalent[_qp]=old_p+dp;_plastic[_qp]=_old_plastic[_qp]+sign*dp;
    _elastic_stretch[_qp]=strain-_plastic[_qp];
    _axial_stress[_qp]=trial-sign*E*dp;
    _e_over_l[_qp]=tangent/length;
  }
}
