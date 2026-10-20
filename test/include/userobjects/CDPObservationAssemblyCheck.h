#pragma once
#include "GeneralUserObject.h"
/** Development-only zero-time assembly exercise; never present in formal inputs. */
class CDPObservationAssemblyCheck : public GeneralUserObject
{
public:
  static InputParameters validParams();
  CDPObservationAssemblyCheck(const InputParameters & p) : GeneralUserObject(p) {}
  void initialize() override {}
  void execute() override;
  void finalize() override {}
};
