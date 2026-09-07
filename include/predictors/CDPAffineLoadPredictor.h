#pragma once
#include "Predictor.h"
class Function;
/** Opt-in initial guess only; no constitutive properties or equilibrium changes. */
class CDPAffineLoadPredictor : public Predictor
{
public:
  static InputParameters validParams();
  CDPAffineLoadPredictor(const InputParameters &);
  bool shouldApply() override;
  void apply(NumericVector<Number> &) override;
private:
  const bool _enabled;
  const std::string _variable;
  const Function & _top_displacement;
  const Real _z_bottom, _z_top;
};
