#pragma once
#include "Exodus.h"
#include "CSV.h"

/** Accepted-step output anchored to ELC marks, without forcing solution times. */
template <typename Base>
class AqueductScheduledOutput : public Base
{
public:
  static InputParameters validParams();
  AqueductScheduledOutput(const InputParameters &);
  void outputStep(const ExecFlagType &) override;
private:
  const Real _origin,_end,_period,_mark_tolerance;
  long long & _last_mark;
};
using AqueductScheduledExodus=AqueductScheduledOutput<Exodus>;
using AqueductScheduledCSV=AqueductScheduledOutput<CSV>;
