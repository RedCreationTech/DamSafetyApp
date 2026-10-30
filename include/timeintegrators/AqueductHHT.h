#pragma once
#include "NewmarkBeta.h"

/** Normalized HHT: M*a/(1+alpha)+F_current-alpha*F_accepted/(1+alpha).
 * Every inertia kernel must use the same hht_alpha. No per-force time weighting.
 */
class AqueductHHT : public NewmarkBeta
{
public:
  static InputParameters validParams();
  AqueductHHT(const InputParameters &);
  void initialSetup() override;
  void init() override;
  void preStep() override;
  void computeTimeDerivatives() override;
  void postResidual(NumericVector<Number> &) override;
  void postStep() override;
  void computeADTimeDerivatives(ADReal &, const dof_id_type &, ADReal &) const override;
  int order() override { return 2; }
protected:
  Real duDotDuCoeff() const override;
  bool quasiStatic() const;
  void clearAlgebraicHistory();
  const bool _quasistatic_initialization;
  Real _hht_alpha, _start_time;
  NumericVector<Number> * _old_forces;
  std::vector<dof_id_type> _algebraic_dofs;
};
