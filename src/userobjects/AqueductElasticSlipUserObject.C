#include "AqueductElasticSlipUserObject.h"
#include "MooseVariableFE.h"
#include "AutomaticMortarGeneration.h"
#include "FEProblemBase.h"

registerMooseObject("DamSafetyApp", AqueductElasticSlipUserObject);

InputParameters AqueductElasticSlipUserObject::validParams()
{
  auto p = LMWeightedVelocitiesUserObject::validParams();
  p.addClassDescription("Three-dimensional mortar elastic-slip history using actual displacement increments");
  return p;
}
AqueductElasticSlipUserObject::AqueductElasticSlipUserObject(const InputParameters & p)
  : WeightedGapUserObject(p), LMWeightedVelocitiesUserObject(p),
    _old_elastic(declareRestartableData<History>("aq_old_elastic")),
    _trial_elastic(declareRestartableData<History>("aq_trial_elastic")),
    _history_step(declareRestartableData<unsigned int>("aq_history_step", 0))
{
  if (!_has_disp_z)
    paramError("disp_z", "This source-model implementation requires three displacement components");
  const MooseVariable * vars[3] = {_disp_x_var, _disp_y_var, _disp_z_var};
  for (unsigned j = 0; j < 3; ++j)
  {
    _secondary_old[j] = &vars[j]->slnOld();
    _primary_old[j] = &vars[j]->slnOldNeighbor();
  }
}
void AqueductElasticSlipUserObject::timestepSetup()
{
  LMWeightedVelocitiesUserObject::timestepSetup();
  // Rejected attempts retain the same step number: never commit their trial history.
  if (_history_step != _fe_problem.timeStep())
  {
    _old_elastic = _trial_elastic;
    _history_step = _fe_problem.timeStep();
  }
}
void AqueductElasticSlipUserObject::initialize()
{
  LMWeightedVelocitiesUserObject::initialize();
  _trial_elastic.clear(); // Nodes without a current projection lose open-contact history.
}
std::array<Real, 2> AqueductElasticSlipUserObject::oldElastic(dof_id_type node) const
{
  const auto i = _old_elastic.find(node);
  return i == _old_elastic.end() ? std::array<Real, 2>{0., 0.} : i->second;
}
void AqueductElasticSlipUserObject::setTrial(dof_id_type node, const std::array<Real, 2> & v) const
{
  _trial_elastic[node] = v;
}
void AqueductElasticSlipUserObject::computeQpProperties()
{
  WeightedGapUserObject::computeQpProperties();
  const auto & pmap = amg().getPrimaryIpToLowerElementMap(
      *_lower_primary_elem, *_lower_primary_elem->interior_parent(), *_lower_secondary_elem);
  const auto & smap = amg().getSecondaryIpToLowerElementMap(*_lower_secondary_elem);
  std::array<const MooseVariable *, 3> vars = {_disp_x_var, _disp_y_var, _disp_z_var};
  std::array<ADReal, 3> primary = {_primary_disp_x[_qp], _primary_disp_y[_qp], (*_primary_disp_z)[_qp]};
  std::array<ADReal, 3> secondary = {_secondary_disp_x[_qp], _secondary_disp_y[_qp], (*_secondary_disp_z)[_qp]};
  trimInteriorNodeDerivatives(pmap, vars, primary, false);
  trimInteriorNodeDerivatives(smap, vars, secondary, true);
  for (unsigned j = 0; j < 3; ++j)
    _qp_real_tangential_velocity_nodal(j) =
        (secondary[j] - (*_secondary_old[j])[_qp] - primary[j] + (*_primary_old[j])[_qp]) / _fe_problem.dt();
  _qp_tangential_velocity_nodal = _qp_real_tangential_velocity_nodal * (_JxW_msm[_qp] * _coord[_qp]);
}
