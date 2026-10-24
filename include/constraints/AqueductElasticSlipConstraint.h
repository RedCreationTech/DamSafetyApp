#pragma once
#include "ComputeFrictionalForceLMMechanicalContact.h"
#include "AqueductElasticSlipUserObject.h"

/** Native hard normal mortar contact plus the source elastic-stick circular return. */
class AqueductElasticSlipConstraint : public ComputeFrictionalForceLMMechanicalContact
{
public:
  static InputParameters validParams();
  AqueductElasticSlipConstraint(const InputParameters & parameters);
protected:
  void enforceConstraintOnDof3d(const DofObject * const dof) override;
  const AqueductElasticSlipUserObject & _elastic_uo;
  const Real _critical_slip;
};
