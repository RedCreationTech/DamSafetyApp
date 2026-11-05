#pragma once
#include "CDPObservationConstantDT.h"
#include <fstream>

/** TEST ONLY: bounded real Transient rejection, not a production CDP failure. */
class DS04GAcceptedIntervalTestStepper : public CDPObservationConstantDT
{
public:
  static InputParameters validParams();
  DS04GAcceptedIntervalTestStepper(const InputParameters &);
  void step() override;
  void acceptStep() override;
  void rejectStep() override;
  void postExecute() override;
private:
  void record(const char *);
  void verifyObservation() const;
  void verifyMaterialObservation() const;
  Real _last_accepted = 0;
  bool _injected = false, _cutback_accepted = false, _after_cutback = false;
  const bool _check_observation;
  const bool _check_material_observation;
  unsigned _accepted = 0;
  std::ofstream _trace;
};
