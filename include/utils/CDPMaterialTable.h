#pragma once

#include "MooseTypes.h"

#include <cstddef>
#include <string>
#include <vector>

/**
 * Validated, piecewise-linear Abaqus CDP material input tables.
 *
 * This class only owns input semantics. It does not implement a yield surface,
 * a flow potential, damage evolution, or stress integration.
 */
class CDPMaterialTable
{
public:
  enum class Branch
  {
    COMPRESSION,
    TENSION
  };

  struct Sample
  {
    Real value;
    Real left_derivative;
    Real right_derivative;
  };

  struct Response
  {
    // Nominal (Cauchy) uniaxial stress from the input table, NOT effective cohesion.
    Sample stress;
    Sample damage;
  };

  CDPMaterialTable(const std::string & compression_hardening_file,
                   const std::string & compression_damage_file,
                   const std::string & tension_stiffening_file,
                   const std::string & tension_damage_file,
                   Real youngs_modulus,
                   Real damage_upper_bound = 1.0 - 1.0e-12);

  Response response(Branch branch, Real strain_measure) const;
  Response responseByEquivalentPlasticStrain(Branch branch, Real equivalent_plastic_strain) const;
  /**
   * Effective uniaxial cohesion sigma_bar = sigma / (1-d) for yield evaluation.
   * Query nominal stress and damage at the same kappa, then apply the quotient
   * rule to both one-sided derivatives. Original input tables are unchanged.
   */
  Sample effectiveStrengthByEquivalentPlasticStrain(Branch branch,
                                                   Real equivalent_plastic_strain) const;
  std::size_t stressPointCount(Branch branch) const;
  std::size_t damagePointCount(Branch branch) const;

  const std::vector<Real> & stressAbscissa(Branch branch) const;
  const std::vector<Real> & stressValues(Branch branch) const;
  const std::vector<Real> & damageAbscissa(Branch branch) const;
  const std::vector<Real> & damageValues(Branch branch) const;
  const std::vector<Real> & equivalentPlasticStrain(Branch branch) const;

private:
  struct Table
  {
    std::string file;
    std::vector<Real> abscissa;
    std::vector<Real> values;
    std::vector<Real> slopes;
  };

  static Table load(const std::string & file,
                    const std::string & value_column,
                    const std::string & abscissa_column,
                    bool damage,
                    Real damage_upper_bound);
  static Sample evaluate(const Table & table, Real abscissa);
  static std::vector<Real> deriveEquivalentPlasticStrain(const Table & stress,
                                                         const Table & damage,
                                                         Real youngs_modulus,
                                                         const std::string & label);
  static Table reparameterizeByEquivalentPlasticStrain(const Table & stress,
                                                       const Table & damage,
                                                       const std::vector<Real> & plastic_strain,
                                                       bool damage_values,
                                                       const std::string & label);

  const Table & stressTable(Branch branch) const;
  const Table & damageTable(Branch branch) const;

  const Table _compression_hardening;
  const Table _compression_damage;
  const Table _tension_stiffening;
  const Table _tension_damage;
  const std::vector<Real> _compression_plastic_strain;
  const std::vector<Real> _tension_plastic_strain;
  const Table _compression_stress_by_plastic_strain;
  const Table _compression_damage_by_plastic_strain;
  const Table _tension_stress_by_plastic_strain;
  const Table _tension_damage_by_plastic_strain;
};
