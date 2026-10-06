#pragma once
#include "TrussMaterial.h"
class SmallStrainPlasticTruss : public TrussMaterial
{
public:
  static InputParameters validParams();
  SmallStrainPlasticTruss(const InputParameters & p);
  void computeProperties() override;
protected:
  void initQpStatefulProperties() override;
  void computeQpStrain() override {}
  void computeQpStress() override {}
  const Real _yield, _ultimate, _knot;
  const MaterialProperty<Real> & _old_strain, & _old_stress;
  MaterialProperty<Real> & _plastic, & _equivalent;
  const MaterialProperty<Real> & _old_plastic, & _old_equivalent;
};
