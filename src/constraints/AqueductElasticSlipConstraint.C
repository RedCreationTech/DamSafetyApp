#include "AqueductElasticSlipConstraint.h"
#include "AqueductElasticSlip.h"
#include "Assembly.h"
#include "MooseVariable.h"

registerMooseObject("DamSafetyApp", AqueductElasticSlipConstraint);

InputParameters AqueductElasticSlipConstraint::validParams()
{
  auto p = ComputeFrictionalForceLMMechanicalContact::validParams();
  p.addRequiredRangeCheckedParam<Real>("critical_slip", "critical_slip>0", "Physical allowable elastic slip in metres, not a search tolerance");
  p.addClassDescription("Source constant-mu elastic stick/slip law with native hard normal mortar contact");
  return p;
}
AqueductElasticSlipConstraint::AqueductElasticSlipConstraint(const InputParameters & p)
  : ComputeFrictionalForceLMMechanicalContact(p),
    _elastic_uo(getUserObject<AqueductElasticSlipUserObject>("weighted_velocities_uo")),
    _critical_slip(getParam<Real>("critical_slip"))
{
  if (!_3d || _has_friction_function || !getParam<bool>("normalize_c"))
    mooseError("AqueductElasticSlipConstraint requires 3D, constant mu and normalize_c=true");
  if (&_weighted_gap_uo != &_weighted_velocities_uo)
    mooseError("The same AqueductElasticSlipUserObject must provide gap and slip");
}
void AqueductElasticSlipConstraint::enforceConstraintOnDof3d(const DofObject * const dof)
{
  ComputeWeightedGapLMMechanicalContact::enforceConstraintOnDof(dof);
  const auto normal_index = dof->dof_number(_sys.number(), _var->number(), 0);
  ADReal pressure = (*_sys.currentSolution())(normal_index);
  Moose::derivInsert(pressure.derivatives(), normal_index, 1.);
  std::array<ADReal, 2> increment;
  for (unsigned j = 0; j < 2; ++j)
    increment[j] = *_tangential_vel_ptr[j] * _dt / *_normalization_ptr;
  const auto state = AqueductFriction::update(_elastic_uo.oldElastic(dof->id()), increment,
                                               pressure, _mu, _critical_slip);
  _elastic_uo.setTrial(dof->id(), {MetaPhysicL::raw_value(state.elastic[0]),
                                  MetaPhysicL::raw_value(state.elastic[1])});
  for (unsigned j = 0; j < 2; ++j)
  {
    const auto index = dof->dof_number(_sys.number(), _friction_vars[j]->number(), 0);
    ADReal traction = (*_sys.currentSolution())(index);
    Moose::derivInsert(traction.derivatives(), index, 1.);
    addResidualsAndJacobian(_assembly, std::array<ADReal, 1>{traction - state.traction[j]},
                            std::array<dof_id_type, 1>{index}, _friction_vars[j]->scalingFactor());
  }
}
