#pragma once
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include <set>

/** Read-only, explicitly gated capture for same-state assembly diagnostics. */
class CDPAssemblyProbe : public Material
{
public:
  static InputParameters validParams();
  CDPAssemblyProbe(const InputParameters & p);
  static void beginCapture(const std::string & path, bool tangent);
  static void endCapture();
protected:
  void computeQpProperties() override;
  const MaterialProperty<RankFourTensor> & _probe_tangent;
  const MaterialProperty<RankTwoTensor> & _probe_stress;
  const MaterialProperty<RankTwoTensor> & _probe_strain;
  const MaterialProperty<Real> & _probe_kappa;
  const MaterialProperty<Real> & _probe_damage;
  const MaterialProperty<Real> & _probe_substeps;
  const MaterialProperty<Real> & _probe_depth;
  const MaterialProperty<Real> & _probe_fallbacks;
  // References only: the observer never updates a constitutive property.
  std::vector<const MaterialProperty<RankTwoTensor> *> _history_tensors;
  std::vector<const MaterialProperty<Real> *> _history_scalars;
  std::set<dof_id_type> _probe_elements;
};
