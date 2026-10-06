#pragma once
#include "ComputeIncrementalStrain.h"
#include "UniformHex8.h"
#include <unordered_map>

// Small-strain uniform-strain HEX8 with one constitutive point. Restricted
// to the undisplaced reference mesh; this does not implement finite strain.
class ComputeUniformHex8IncrementalStrain : public ComputeIncrementalStrain
{
public:
  static InputParameters validParams();
  ComputeUniformHex8IncrementalStrain(const InputParameters & parameters);
protected:
  void computeTotalStrainIncrement(RankTwoTensor & increment) override;
  std::array<MooseVariable *, 3> _nodal_displacements;
  std::unordered_map<dof_id_type, UniformHex8::Operators> _operators;
};
