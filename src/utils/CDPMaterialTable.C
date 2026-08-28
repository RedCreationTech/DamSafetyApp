#include "CDPMaterialTable.h"

#include "DelimitedFileReader.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace
{
[[noreturn]] void
tableError(const std::string & message)
{
  throw std::invalid_argument("CDPMaterialTable: " + message);
}
}

CDPMaterialTable::CDPMaterialTable(const std::string & compression_hardening_file,
                                   const std::string & compression_damage_file,
                                   const std::string & tension_stiffening_file,
                                   const std::string & tension_damage_file,
                                   const Real youngs_modulus,
                                   const Real damage_upper_bound)
  : _compression_hardening(load(
        compression_hardening_file, "stress_pa", "inelastic_strain", false, damage_upper_bound)),
    _compression_damage(
        load(compression_damage_file, "damage_c", "inelastic_strain", true, damage_upper_bound)),
    _tension_stiffening(
        load(tension_stiffening_file, "stress_pa", "cracking_strain", false, damage_upper_bound)),
    _tension_damage(
        load(tension_damage_file, "damage_t", "cracking_strain", true, damage_upper_bound)),
    _compression_plastic_strain(deriveEquivalentPlasticStrain(
        _compression_hardening, _compression_damage, youngs_modulus, "compression")),
    _tension_plastic_strain(deriveEquivalentPlasticStrain(
        _tension_stiffening, _tension_damage, youngs_modulus, "tension")),
    _compression_stress_by_plastic_strain(reparameterizeByEquivalentPlasticStrain(
        _compression_hardening, _compression_damage, _compression_plastic_strain, false, "compression")),
    _compression_damage_by_plastic_strain(reparameterizeByEquivalentPlasticStrain(
        _compression_hardening, _compression_damage, _compression_plastic_strain, true, "compression")),
    _tension_stress_by_plastic_strain(reparameterizeByEquivalentPlasticStrain(
        _tension_stiffening, _tension_damage, _tension_plastic_strain, false, "tension")),
    _tension_damage_by_plastic_strain(reparameterizeByEquivalentPlasticStrain(
        _tension_stiffening, _tension_damage, _tension_plastic_strain, true, "tension"))
{
  if (!std::isfinite(youngs_modulus) || youngs_modulus <= 0.0)
    tableError("Young's modulus must be finite and positive");
  if (!std::isfinite(damage_upper_bound) || damage_upper_bound <= 0.0 || damage_upper_bound >= 1.0)
    tableError("damage_upper_bound must be finite and in (0, 1)");
}

CDPMaterialTable::Table
CDPMaterialTable::load(const std::string & file,
                       const std::string & value_column,
                       const std::string & abscissa_column,
                       const bool damage,
                       const Real damage_upper_bound)
{
  MooseUtils::DelimitedFileReader reader(file);
  reader.setHeaderFlag(MooseUtils::DelimitedFileReader::HeaderFlag::ON);
  reader.setFormatFlag(MooseUtils::DelimitedFileReader::FormatFlag::COLUMNS);
  reader.read();

  if (reader.getData().size() != 2)
    tableError(file + " must contain exactly two columns");

  const auto & names = reader.getNames();
  if (std::find(names.begin(), names.end(), value_column) == names.end() ||
      std::find(names.begin(), names.end(), abscissa_column) == names.end())
    tableError(file + " must contain columns '" + value_column + "' and '" + abscissa_column + "'");

  Table table;
  table.file = file;
  table.values = reader.getData(value_column);
  table.abscissa = reader.getData(abscissa_column);

  if (table.values.empty() || table.values.size() != table.abscissa.size())
    tableError(file + " must contain at least one complete data row");

  for (std::size_t i = 0; i < table.values.size(); ++i)
  {
    const Real x = table.abscissa[i];
    const Real y = table.values[i];
    if (!std::isfinite(x) || !std::isfinite(y))
      tableError(file + " contains a non-finite value at data row " + std::to_string(i + 1));
    if (x < 0.0)
      tableError(file + " contains a negative strain abscissa at data row " +
                 std::to_string(i + 1));
    if (i > 0 && x <= table.abscissa[i - 1])
      tableError(file + " strain abscissa must be strictly increasing at data row " +
                 std::to_string(i + 1));
    if (damage)
    {
      if (y < 0.0 || y > damage_upper_bound)
        tableError(file + " damage must be in [0, damage_upper_bound] at data row " +
                   std::to_string(i + 1));
    }
    else if (y <= 0.0)
      tableError(file + " stress magnitude must be positive at data row " + std::to_string(i + 1));
  }

  if (table.abscissa.front() != 0.0)
    tableError(file + " first strain abscissa must be zero");
  if (damage && table.values.front() != 0.0)
    tableError(file + " first damage value must be zero");

  table.slopes.reserve(table.values.size() > 1 ? table.values.size() - 1 : 0);
  for (std::size_t i = 1; i < table.values.size(); ++i)
    table.slopes.push_back((table.values[i] - table.values[i - 1]) /
                           (table.abscissa[i] - table.abscissa[i - 1]));

  return table;
}

CDPMaterialTable::Sample
CDPMaterialTable::evaluate(const Table & table, const Real abscissa)
{
  if (!std::isfinite(abscissa))
    tableError("query abscissa must be finite");

  if (table.values.size() == 1)
    return {table.values.front(), 0.0, 0.0};

  if (abscissa < table.abscissa.front())
    return {table.values.front(), 0.0, 0.0};
  if (abscissa == table.abscissa.front())
    return {table.values.front(), 0.0, table.slopes.front()};
  if (abscissa > table.abscissa.back())
    return {table.values.back(), 0.0, 0.0};
  if (abscissa == table.abscissa.back())
    return {table.values.back(), table.slopes.back(), 0.0};

  const auto upper = std::lower_bound(table.abscissa.begin(), table.abscissa.end(), abscissa);
  const std::size_t upper_index = std::distance(table.abscissa.begin(), upper);
  if (*upper == abscissa)
    return {table.values[upper_index], table.slopes[upper_index - 1], table.slopes[upper_index]};

  const std::size_t lower_index = upper_index - 1;
  const Real slope = table.slopes[lower_index];
  const Real value = table.values[lower_index] + slope * (abscissa - table.abscissa[lower_index]);
  return {value, slope, slope};
}

std::vector<Real>
CDPMaterialTable::deriveEquivalentPlasticStrain(const Table & stress,
                                                const Table & damage,
                                                const Real youngs_modulus,
                                                const std::string & label)
{
  if (!std::isfinite(youngs_modulus) || youngs_modulus <= 0.0)
    tableError("Young's modulus must be finite and positive");

  std::vector<Real> result;
  result.reserve(stress.values.size());
  Real previous = 0.0;
  for (std::size_t i = 0; i < stress.values.size(); ++i)
  {
    const Real strain_measure = stress.abscissa[i];
    const Real sigma = stress.values[i];
    const Real d = evaluate(damage, strain_measure).value;
    const Real plastic_strain = strain_measure - d / (1.0 - d) * sigma / youngs_modulus;
    const Real tolerance = 1.0e-14 * std::max(1.0, std::abs(strain_measure));

    if (plastic_strain < -tolerance)
    {
      std::ostringstream message;
      message << label << " equivalent plastic strain is negative at data row " << i + 1 << ": "
              << plastic_strain;
      tableError(message.str());
    }

    const Real clamped = std::max(0.0, plastic_strain);
    if (!result.empty() && clamped + tolerance < previous)
    {
      std::ostringstream message;
      message << label << " equivalent plastic strain decreases at data row " << i + 1 << ": "
              << clamped << " < " << previous;
      tableError(message.str());
    }
    result.push_back(clamped);
    previous = clamped;
  }
  return result;
}

CDPMaterialTable::Table
CDPMaterialTable::reparameterizeByEquivalentPlasticStrain(
    const Table & stress,
    const Table & damage,
    const std::vector<Real> & plastic_strain,
    const bool damage_values,
    const std::string & label)
{
  if (plastic_strain.size() != stress.values.size())
    tableError(label + " plastic strain and stress table sizes differ");

  Table result;
  result.file = label + (damage_values ? " damage" : " stress") +
                " by equivalent plastic strain";
  result.abscissa = plastic_strain;
  result.values.reserve(stress.values.size());
  for (std::size_t i = 0; i < stress.values.size(); ++i)
  {
    if (i > 0 && result.abscissa[i] <= result.abscissa[i - 1])
      tableError(label + " equivalent plastic strain must be strictly increasing at data row " +
                 std::to_string(i + 1));
    result.values.push_back(damage_values ? evaluate(damage, stress.abscissa[i]).value
                                          : stress.values[i]);
  }

  result.slopes.reserve(result.values.size() > 1 ? result.values.size() - 1 : 0);
  for (std::size_t i = 1; i < result.values.size(); ++i)
    result.slopes.push_back((result.values[i] - result.values[i - 1]) /
                            (result.abscissa[i] - result.abscissa[i - 1]));
  return result;
}

const CDPMaterialTable::Table &
CDPMaterialTable::stressTable(const Branch branch) const
{
  return branch == Branch::COMPRESSION ? _compression_hardening : _tension_stiffening;
}

const CDPMaterialTable::Table &
CDPMaterialTable::damageTable(const Branch branch) const
{
  return branch == Branch::COMPRESSION ? _compression_damage : _tension_damage;
}

CDPMaterialTable::Response
CDPMaterialTable::response(const Branch branch, const Real strain_measure) const
{
  return {evaluate(stressTable(branch), strain_measure),
          evaluate(damageTable(branch), strain_measure)};
}

CDPMaterialTable::Response
CDPMaterialTable::responseByEquivalentPlasticStrain(const Branch branch,
                                                    const Real equivalent_plastic_strain) const
{
  if (branch == Branch::COMPRESSION)
    return {evaluate(_compression_stress_by_plastic_strain, equivalent_plastic_strain),
            evaluate(_compression_damage_by_plastic_strain, equivalent_plastic_strain)};
  return {evaluate(_tension_stress_by_plastic_strain, equivalent_plastic_strain),
          evaluate(_tension_damage_by_plastic_strain, equivalent_plastic_strain)};
}

std::size_t
CDPMaterialTable::stressPointCount(const Branch branch) const
{
  return stressTable(branch).values.size();
}

std::size_t
CDPMaterialTable::damagePointCount(const Branch branch) const
{
  return damageTable(branch).values.size();
}

const std::vector<Real> &
CDPMaterialTable::stressAbscissa(const Branch branch) const
{
  return stressTable(branch).abscissa;
}

const std::vector<Real> &
CDPMaterialTable::stressValues(const Branch branch) const
{
  return stressTable(branch).values;
}

const std::vector<Real> &
CDPMaterialTable::damageAbscissa(const Branch branch) const
{
  return damageTable(branch).abscissa;
}

const std::vector<Real> &
CDPMaterialTable::damageValues(const Branch branch) const
{
  return damageTable(branch).values;
}

const std::vector<Real> &
CDPMaterialTable::equivalentPlasticStrain(const Branch branch) const
{
  return branch == Branch::COMPRESSION ? _compression_plastic_strain : _tension_plastic_strain;
}
