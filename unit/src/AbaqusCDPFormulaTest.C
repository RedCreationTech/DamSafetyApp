#include "AbaqusCDPFormula.h"

#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <stdexcept>

using AbaqusCDPFormula::SymmetricTensor;

TEST(AbaqusCDPFormula, ComputesPressureMisesPrincipalStressAndTensionWeight)
{
  const auto hydrostatic =
      AbaqusCDPFormula::stressInvariants({-30.0, -30.0, -30.0, 0.0, 0.0, 0.0});
  EXPECT_DOUBLE_EQ(hydrostatic.pressure, 30.0);
  EXPECT_DOUBLE_EQ(hydrostatic.mises, 0.0);
  EXPECT_DOUBLE_EQ(hydrostatic.tension_weight, 0.0);

  const auto tension = AbaqusCDPFormula::stressInvariants({10.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  EXPECT_NEAR(tension.pressure, -10.0 / 3.0, 1.0e-14);
  EXPECT_NEAR(tension.mises, 10.0, 1.0e-14);
  EXPECT_DOUBLE_EQ(tension.maximum_principal_stress, 10.0);
  EXPECT_DOUBLE_EQ(tension.tension_weight, 1.0);

  const auto mixed = AbaqusCDPFormula::stressInvariants({-10.0, 5.0, 0.0, 0.0, 0.0, 0.0});
  EXPECT_NEAR(mixed.tension_weight, 1.0 / 3.0, 1.0e-14);

  const auto sheared = AbaqusCDPFormula::stressInvariants({4.0, -2.0, 1.0, 3.0, -1.0, 2.0});
  EXPECT_LE(sheared.principal_stress[0], sheared.principal_stress[1]);
  EXPECT_LE(sheared.principal_stress[1], sheared.principal_stress[2]);
  EXPECT_NEAR(sheared.principal_stress[0] + sheared.principal_stress[1] +
                  sheared.principal_stress[2],
              3.0,
              1.0e-13);
}

TEST(AbaqusCDPFormula, ReproducesUniaxialAndBiaxialYieldCalibrationPoints)
{
  constexpr double compression_strength = 30.0;
  constexpr double tension_strength = 3.0;
  constexpr double biaxial_ratio = 1.16;
  constexpr double kc = 0.667;

  EXPECT_NEAR(AbaqusCDPFormula::yieldFunction({-compression_strength, 0, 0, 0, 0, 0},
                                              compression_strength,
                                              tension_strength,
                                              biaxial_ratio,
                                              kc),
              0.0,
              1.0e-12);
  EXPECT_NEAR(AbaqusCDPFormula::yieldFunction({tension_strength, 0, 0, 0, 0, 0},
                                              compression_strength,
                                              tension_strength,
                                              biaxial_ratio,
                                              kc),
              0.0,
              1.0e-12);
  const double biaxial_strength = biaxial_ratio * compression_strength;
  EXPECT_NEAR(AbaqusCDPFormula::yieldFunction({-biaxial_strength,
                                               -biaxial_strength,
                                               0,
                                               0,
                                               0,
                                               0},
                                              compression_strength,
                                              tension_strength,
                                              biaxial_ratio,
                                              kc),
              0.0,
              1.0e-12);
}

TEST(AbaqusCDPFormula, FlowPotentialGradientMatchesFiniteDifferences)
{
  const SymmetricTensor stress = {4.0e6, -2.0e6, 1.0e6, 0.5e6, -0.25e6, 0.75e6};
  const auto evaluated = AbaqusCDPFormula::flowPotential(stress, 36.31, 0.1, 2.54981e6);
  constexpr double perturbation = 1.0;
  for (std::size_t i = 0; i < stress.size(); ++i)
  {
    auto plus = stress;
    auto minus = stress;
    plus[i] += perturbation;
    minus[i] -= perturbation;
    const double finite_difference =
        (AbaqusCDPFormula::flowPotential(plus, 36.31, 0.1, 2.54981e6).value -
         AbaqusCDPFormula::flowPotential(minus, 36.31, 0.1, 2.54981e6).value) /
        (2.0 * perturbation);
    EXPECT_NEAR(evaluated.gradient[i], finite_difference, 2.0e-9);
  }
}

TEST(AbaqusCDPFormula, CombinesDamageWithIndependentRecoveryFactors)
{
  const auto pure_tension = AbaqusCDPFormula::combineDamage(0.8, 0.2, 1.0, 0.0, 1.0);
  EXPECT_DOUBLE_EQ(pure_tension.tensile_recovery_factor, 0.0);
  EXPECT_DOUBLE_EQ(pure_tension.damage, 0.2);

  const auto pure_compression = AbaqusCDPFormula::combineDamage(0.8, 0.2, 0.0, 1.0, 0.0);
  EXPECT_DOUBLE_EQ(pure_compression.compressive_recovery_factor, 0.0);
  EXPECT_DOUBLE_EQ(pure_compression.damage, 0.8);

  const auto mixed = AbaqusCDPFormula::combineDamage(0.4, 0.25, 0.8, 0.6, 0.25);
  EXPECT_GE(mixed.damage, 0.0);
  EXPECT_LE(mixed.damage, 1.0);
  EXPECT_NEAR(mixed.damage, 1.0 - (1.0 - 0.8 * 0.4) * (1.0 - 0.55 * 0.25), 1.0e-14);
}

TEST(AbaqusCDPFormula, AppliesDuvautLionsBackwardEulerUpdate)
{
  EXPECT_DOUBLE_EQ(AbaqusCDPFormula::duvautLionsUpdate(0.0, 1.0, 0.5, 0.5), 0.5);
  EXPECT_DOUBLE_EQ(AbaqusCDPFormula::duvautLionsUpdate(0.2, 1.0, 0.0, 0.5), 0.2);
  EXPECT_DOUBLE_EQ(AbaqusCDPFormula::duvautLionsUpdate(0.2, 1.0, 0.5, 0.0), 1.0);

  const SymmetricTensor old_state = {0, 1, 2, 3, 4, 5};
  const SymmetricTensor backbone = {2, 3, 4, 5, 6, 7};
  const auto updated = AbaqusCDPFormula::duvautLionsUpdate(old_state, backbone, 1.0, 1.0);
  for (std::size_t i = 0; i < updated.size(); ++i)
    EXPECT_DOUBLE_EQ(updated[i], old_state[i] + 1.0);
}

TEST(AbaqusCDPFormula, RejectsInvalidOfficialParameters)
{
  EXPECT_THROW(AbaqusCDPFormula::yieldCoefficients(1.16, 0.5, 30.0, 3.0),
               std::invalid_argument);
  EXPECT_THROW(AbaqusCDPFormula::yieldCoefficients(1.16, 1.01, 30.0, 3.0),
               std::invalid_argument);
  EXPECT_THROW(AbaqusCDPFormula::flowPotential({0, 0, 0, 0, 0, 0}, 90.0, 0.1, 3.0),
               std::invalid_argument);
  EXPECT_THROW(AbaqusCDPFormula::duvautLionsUpdate(0.0, 1.0, -1.0, 0.5),
               std::invalid_argument);
}
