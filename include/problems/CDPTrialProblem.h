#pragma once
#include "FEProblem.h"
#include <petscsnes.h>
/** Passive observation of real residual evaluations and SNES accepted iterates. */
class CDPTrialProblem : public FEProblem
{
public:
  static InputParameters validParams();
  CDPTrialProblem(const InputParameters & p);
  void computeResidual(const NumericVector<Number> &, NumericVector<Number> &, unsigned int) override;
  static PetscErrorCode monitor(SNES, PetscInt, PetscReal, void *);
private:
  bool capturing() const;
  void auditDirections(SNES);
  bool _directions_done = false;
  unsigned int _evaluation = 0;
  const Real _start, _end;
  const std::string _prefix;
};
