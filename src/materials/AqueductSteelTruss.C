#include "AqueductSteelTruss.h"
#include "FEProblem.h"
#include "NonlinearSystemBase.h"
#include "MooseVariable.h"
#include "libmesh/elem.h"

registerMooseObject("DamSafetyApp", AqueductSteelTruss);

InputParameters AqueductSteelTruss::validParams()
{
  auto p=TrussMaterial::validParams();
  p.addRequiredRangeCheckedParam<Real>("yield_stress","yield_stress>0","Initial source yield stress");
  p.addRequiredRangeCheckedParam<Real>("ultimate_stress","ultimate_stress>0","Last source stress; constant extrapolation");
  p.addRequiredRangeCheckedParam<Real>("plastic_strain_at_ultimate","plastic_strain_at_ultimate>0","Last source plastic strain");
  p.set<bool>("use_displaced_mesh")=false;
  return p;
}
AqueductSteelTruss::AqueductSteelTruss(const InputParameters & p)
 :TrussMaterial(p),_fy(getParam<Real>("yield_stress")),_fu(getParam<Real>("ultimate_stress")),
 _limit(getParam<Real>("plastic_strain_at_ultimate")),
 _plastic(declareProperty<Real>(_base_name+"signed_plastic_stretch")),
 _plastic_old(getMaterialPropertyOld<Real>(_base_name+"signed_plastic_stretch")),
 _equivalent(declareProperty<Real>(_base_name+"equivalent_plastic_stretch")),
 _equivalent_old(getMaterialPropertyOld<Real>(_base_name+"equivalent_plastic_stretch"))
{
  if(_ndisp!=3)paramError("displacements","AQ1023 T3D2 requires three displacement variables");
  if(_fu<_fy)paramError("ultimate_stress","Last stress must not be below initial yield");
  if(getParam<bool>("use_displaced_mesh"))paramError("use_displaced_mesh","Source NLGEOM=NO requires the reference mesh");
}
void AqueductSteelTruss::initQpStatefulProperties()
{
  TrussMaterial::initQpStatefulProperties();_plastic[_qp]=0;_equivalent[_qp]=0;
}
void AqueductSteelTruss::computeQpStrain()
{
  auto & system=_fe_problem.getNonlinearSystemBase(0);
  const auto & solution=*system.currentSolution();
  const auto & a=_current_elem->node_ref(0);const auto & b=_current_elem->node_ref(1);
  _total_stretch[_qp]=0;
  for(unsigned j=0;j<_ndisp;++j)
    _total_stretch[_qp]+=(b(j)-a(j))*
      (solution(b.dof_number(system.number(),_disp_var[j]->number(),0))-
       solution(a.dof_number(system.number(),_disp_var[j]->number(),0)))/(_origin_length*_origin_length);
}
void AqueductSteelTruss::computeQpStress()
{
  auto s=AqueductTruss::update(_total_stretch[_qp],_plastic_old[_qp],_equivalent_old[_qp],
                              _youngs_modulus[_qp],_fy,_fu,_limit);
  _axial_stress[_qp]=s.stress;_plastic[_qp]=s.plastic;_equivalent[_qp]=s.equivalent;
  _elastic_stretch[_qp]=_total_stretch[_qp]-s.plastic;
  _e_over_l[_qp]=s.tangent/_origin_length;
}
