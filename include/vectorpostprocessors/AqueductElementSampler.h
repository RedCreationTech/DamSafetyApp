#pragma once
#include "ElementValueSampler.h"
#include <set>

/** Direct elemental values selected by immutable libMesh element IDs. */
class AqueductElementSampler : public ElementValueSampler
{
public:
  static InputParameters validParams();
  AqueductElementSampler(const InputParameters &);
  void initialSetup() override;
  void execute() override;
private:
  const std::set<dof_id_type> _ids;
};
