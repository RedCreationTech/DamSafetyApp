#pragma once
#include "FileOutput.h"

/** Opt-in accepted-state diagnostic. Does not alter the CDP constitutive equations. */
class CDPAcceptedStateOutput : public FileOutput
{
public:
  static InputParameters validParams();
  CDPAcceptedStateOutput(const InputParameters & parameters);
  std::string filename() override;
protected:
  void output() override;
  void writeHistory(const std::string & stage);
};
