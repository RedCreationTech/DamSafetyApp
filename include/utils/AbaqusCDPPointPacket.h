#pragma once

#include "AbaqusCDPSubstepIntegrator.h"

#include <cstdint>
#include <istream>
#include <limits>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

/** Self-contained diagnostic input. Floating values use IEEE-754 bits, not rounded logs. */
struct AbaqusCDPPointPacket
{
  std::uint64_t element = 0, step = 0;
  unsigned int qp = 0, rank = 0, thread = 0;
  unsigned int side = std::numeric_limits<unsigned int>::max();
  double time = 0., dt = 0.;
  std::array<double, 3> position = {};
  std::string material, expected_error;
  AbaqusCDPLocalIntegrator::Parameters local = {};
  AbaqusCDPStateIntegrator::Parameters state = {};
  AbaqusCDPSubstepIntegrator::Parameters substeps;
  AbaqusCDPSubstepIntegrator::SymmetricTensor old_strain = {}, new_strain = {};
  AbaqusCDPSubstepIntegrator::State old;
  // Each pair is {original CSV value, original strain abscissa}.
  std::array<std::vector<std::pair<double, double>>, 4> tables;

  void setTables(const CDPMaterialTable & table);
  void write(std::ostream & output) const;
  static AbaqusCDPPointPacket read(std::istream & input);
  std::array<std::string, 4> writeTables(const std::string & directory) const;
};
