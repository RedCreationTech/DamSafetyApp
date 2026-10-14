// Developer-only replay of H43 checkpoints, not a finite-element solve.
#include "AbaqusCDPSubstepIntegrator.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

using State = AbaqusCDPStateIntegrator::State;
using Tensor = AbaqusCDPStateIntegrator::SymmetricTensor;

std::array<double, 17> flatten(const State & state)
{
  std::array<double, 17> a{};
  std::copy(state.backbone.plastic_strain.begin(), state.backbone.plastic_strain.end(), a.begin());
  a[6] = state.backbone.tensile_equivalent_plastic_strain;
  a[7] = state.backbone.compressive_equivalent_plastic_strain;
  std::copy(state.viscous_plastic_strain.begin(), state.viscous_plastic_strain.end(), a.begin() + 8);
  a[14] = state.viscous_tension_damage;
  a[15] = state.viscous_compression_damage;
  a[16] = state.viscous_combined_damage;
  return a;
}

int main(int argc, char ** argv)
{
  try
  {
    if (argc != 3)
      throw std::runtime_error("Usage: replay_dam2b_checkpoint material_directory snapshots.csv");
    const std::string p = std::string(argv[1]) + "/";
    constexpr double E = 3.04e10;
    CDPMaterialTable table(p + "compression_hardening.csv", p + "compression_damage.csv",
                          p + "tension_stiffening.csv", p + "tension_damage.csv", E);
    const AbaqusCDPLocalIntegrator local(table, {E, .2, 36.31, .1, 1.16, .667,
                                               40, 1e-9, 1e-7, 1e-6, true, false});
    std::ifstream input(argv[2]);
    if (!input)
      throw std::runtime_error("Cannot open checkpoint snapshots");
    std::string line;
    std::getline(input, line); // schema header is independently checked by the snapshot producer
    std::cout << std::setprecision(17)
              << "abaqus_element,mode,dt,h,stress_replay_absmax_Pa,state_replay_absmax,"
                 "C3312_Pa,C1233_Pa,tangent_column,FD_absolute_error_Pa,FD_relative_error,"
                 "accepted_substeps,old_state_unchanged,selected_analytic_Pa,selected_FD_center_Pa,"
                 "selected_FD_plus_Pa,selected_FD_minus_Pa,base_S33_Pa,plus_substeps,minus_substeps\n";
    std::size_t rows = 0;
    while (std::getline(input, line))
    {
      if (!line.empty() && line.back() == '\r')
        line.pop_back();
      std::istringstream fields(line);
      std::string field;
      std::vector<double> x;
      while (std::getline(fields, field, ','))
      {
        std::size_t consumed;
        const double value = std::stod(field, &consumed);
        if (consumed != field.size() || !std::isfinite(value))
          throw std::runtime_error("Nonfinite or malformed snapshot field");
        x.push_back(value);
      }
      if (x.size() != 56 || !(x[3] > 0) || (x[2] != 0 && x[2] != 1))
        throw std::runtime_error("Invalid snapshot schema");
      Tensor eold{}, enew{}, expected_stress{};
      State old;
      std::copy_n(x.begin() + 4, 6, eold.begin());
      std::copy_n(x.begin() + 10, 6, enew.begin());
      std::copy_n(x.begin() + 16, 6, old.backbone.plastic_strain.begin());
      old.backbone.tensile_equivalent_plastic_strain = x[22];
      old.backbone.compressive_equivalent_plastic_strain = x[23];
      std::copy_n(x.begin() + 24, 6, old.viscous_plastic_strain.begin());
      old.viscous_tension_damage = x[30];
      old.viscous_compression_damage = x[31];
      old.viscous_combined_damage = x[32];
      std::copy_n(x.begin() + 33, 6, expected_stress.begin());
      const auto saved_old = flatten(old);
      const AbaqusCDPStateIntegrator core(local, {1., 1., .0005, 1e-12, x[2] == 1});
      const AbaqusCDPSubstepIntegrator steps(core, {256, 2.5e-5, 1e-8});
      const auto result = steps.integrateLinearized(eold, enew, x[3], old);
      const auto state = flatten(result.result.final_result.state);
      double stress_error = 0, state_error = 0;
      for (unsigned int i = 0; i < 6; ++i)
        stress_error = std::max(stress_error,
                               std::abs(result.result.final_result.cauchy_stress[i] - expected_stress[i]));
      for (unsigned int i = 0; i < 17; ++i)
        state_error = std::max(state_error, std::abs(state[i] - x[39 + i]));
      // Output replication is a prerequisite; do not certify a tangent for the wrong time layer.
      if (stress_error > 1e-4 || state_error > 1e-11 || flatten(old) != saved_old)
        throw std::runtime_error("Checkpoint replay/history identity failed; inspect output before interpreting tangent");
      for (unsigned int col : {2u, 3u})
        for (double h : {1e-8, 1e-9, 1e-10, 1e-12, 1e-14, 1e-16})
        {
          Tensor direction{};
          direction[col] = 1;
          const auto fd = steps.directionalDerivative(eold, enew, x[3], old, direction, h);
          auto plus_strain = enew, minus_strain = enew;
          plus_strain[col] += h;
          minus_strain[col] -= h;
          const auto plus = steps.integrate(eold, plus_strain, x[3], old);
          const auto minus = steps.integrate(eold, minus_strain, x[3], old);
          const unsigned int selected_row = col == 2 ? 3 : 2;
          const double selected_factor = col == 2 ? 1. : .5;
          const double selected_analytic = selected_factor * result.algorithmic_tangent[col][selected_row];
          const double selected_plus = selected_factor * (plus.final_result.cauchy_stress[selected_row] - result.result.final_result.cauchy_stress[selected_row]) / h;
          const double selected_minus = selected_factor * (result.result.final_result.cauchy_stress[selected_row] - minus.final_result.cauchy_stress[selected_row]) / h;
          double error = 0, scale = 0;
          for (unsigned int i = 0; i < 6; ++i)
          {
            error = std::max(error, std::abs(fd[i] - result.algorithmic_tangent[col][i]));
            scale = std::max(scale, std::abs(result.algorithmic_tangent[col][i]));
          }
          std::cout << x[0] << ',' << x[2] << ',' << x[3] << ',' << h << ','
                    << stress_error << ',' << state_error << ','
                    << .5 * result.algorithmic_tangent[3][2] << ','
                    << result.algorithmic_tangent[2][3] << ',' << col << ','
                    << error << ',' << error / std::max(scale, 1.) << ','
                    << result.result.accepted_substeps << ',' << (flatten(old) == saved_old) << ','
                    << selected_analytic << ',' << selected_factor * fd[selected_row] << ','
                    << selected_plus << ',' << selected_minus << ','
                    << result.result.final_result.cauchy_stress[2] << ','
                    << plus.accepted_substeps << ',' << minus.accepted_substeps << '\n';
        }
      ++rows;
    }
    if (!rows)
      throw std::runtime_error("No checkpoint snapshots");
    return 0;
  }
  catch (const std::exception & error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
