#pragma once
#include "LMWeightedVelocitiesUserObject.h"
#include <array>
#include <map>

/** Mortar displacement increments and accepted elastic-slip history, independent of Newmark v*dt. */
class AqueductElasticSlipUserObject : public LMWeightedVelocitiesUserObject
{
public:
  static InputParameters validParams();
  AqueductElasticSlipUserObject(const InputParameters & parameters);
  void timestepSetup() override;
  void initialize() override;
  std::array<Real, 2> oldElastic(dof_id_type node) const;
  void setTrial(dof_id_type node, const std::array<Real, 2> & value) const;
protected:
  void computeQpProperties() override;
  std::array<const VariableValue *, 3> _secondary_old;
  std::array<const VariableValue *, 3> _primary_old;
  using History = std::map<dof_id_type, std::array<Real, 2>>;
  History & _old_elastic;
  History & _trial_elastic;
  unsigned int & _history_step;
};
