#include "CDPObservationMaterial.h"
#include "CDPObservation.h"
#include "FEProblemBase.h"
#include "libmesh/elem.h"
registerMooseObject("DamSafetyApp", CDPObservationMaterial);
InputParameters CDPObservationMaterial::validParams() { return Material::validParams(); }
CDPObservationMaterial::CDPObservationMaterial(const InputParameters & p)
  : Material(p), _tangent(getMaterialProperty<RankFourTensor>("Jacobian_mult")),
    _stress(getMaterialProperty<RankTwoTensor>("stress")),
    _strain(getMaterialProperty<RankTwoTensor>("mechanical_strain")) {}
void CDPObservationMaterial::computeQpProperties()
{
  CDPObservation::PointScope point(_current_elem->id(),_qp,false);
  if(!CDPObservation::active()) return;
  CDPObservation::values("probe-map","qp","x_y_z_JxW_coord",{_q_point[_qp](0),_q_point[_qp](1),_q_point[_qp](2),_JxW[_qp],_coord[_qp]});
  std::array<double,9> stress,strain;
  for(unsigned int c=0;c<9;++c) { stress[c]=_stress[_qp](c/3,c%3); strain[c]=_strain[_qp](c/3,c%3); }
  CDPObservation::array("qp-updates","assembled","stress_row_major",stress);
  CDPObservation::array("qp-updates","assembled","mechanical_strain_row_major",strain);
  if(_fe_problem.currentlyComputingJacobian())
  {
    std::array<double,81> tangent;
    for(unsigned int c=0;c<81;++c) tangent[c]=_tangent[_qp](c/27,(c/9)%3,(c/3)%3,c%3);
    CDPObservation::array("algorithmic-tangent","assembly","RankFour_ijkl",tangent);
  }
}
