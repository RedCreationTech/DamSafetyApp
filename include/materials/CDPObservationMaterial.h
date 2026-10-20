#pragma once
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
/** Reads current assembly quantities only; requests no Old/Older properties. */
class CDPObservationMaterial : public Material
{
public:
  static InputParameters validParams();
  CDPObservationMaterial(const InputParameters &);
protected:
  void computeQpProperties() override;
  const MaterialProperty<RankFourTensor> & _tangent;
  const MaterialProperty<RankTwoTensor> & _stress, & _strain;
};
