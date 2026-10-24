#pragma once

#include "Material.h"
#include "ADKernel.h"
#include "ADTimeKernel.h"
#include "AqueductShellSection.h"

// AQ1023 flat rectangular six-DOF MITC4 engineering surrogate; SMALL strain only.
class AqueductFlatShellStress : public Material
{
public:
  static InputParameters validParams();
  AqueductFlatShellStress(const InputParameters &);
protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;
  std::array<const ADVariableValue *,6> _dofs;
  Real _h,_offset,_E,_nu,_fy,_fu,_limit,_shear_factor,_drill_factor;
  MaterialProperty<std::vector<Real>> & _state;
  const MaterialProperty<std::vector<Real>> & _old_state;
  MaterialProperty<std::vector<Real>> & _strain;
  const MaterialProperty<std::vector<Real>> & _old_strain;
  ADMaterialProperty<RealVectorValue> & _membrane;
  ADMaterialProperty<RealVectorValue> & _moment;
  ADMaterialProperty<RealVectorValue> & _shear;
  ADMaterialProperty<Real> & _drill;
};

class AqueductFlatShell : public ADKernel
{
public:
  static InputParameters validParams();
  AqueductFlatShell(const InputParameters &);
protected:
  ADReal computeQpResidual() override;
  unsigned int _component;
  const ADMaterialProperty<RealVectorValue> & _membrane;
  const ADMaterialProperty<RealVectorValue> & _moment;
  const ADMaterialProperty<RealVectorValue> & _shear;
  const ADMaterialProperty<Real> & _drill;
};

class AqueductFlatShellInertia : public ADTimeKernel
{
public:
  static InputParameters validParams();
  AqueductFlatShellInertia(const InputParameters &);
protected:
  ADReal computeQpResidual() override;
  unsigned int _component;
  Real _rho;
  Real _hht_mass_scale, _hht_start_time;
  AqueductShell::Moments _mass;
  std::array<const ADVariableValue *,6> _accelerations;
};
