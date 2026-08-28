#include "CDPMaterialTable.h"

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace
{
const std::string data_dir = "test/tests/cdp_material_table/data/";

CDPMaterialTable
referenceTable()
{
  return CDPMaterialTable(data_dir + "compression_hardening.csv",
                          data_dir + "compression_damage.csv",
                          data_dir + "tension_stiffening.csv",
                          data_dir + "tension_damage.csv",
                          3.04e10);
}

std::vector<std::pair<Real, Real>>
readPythonBundleTable(const std::string & filename)
{
  std::ifstream stream(filename);
  EXPECT_TRUE(stream.good());
  std::string line;
  std::getline(stream, line);

  std::vector<std::pair<Real, Real>> rows;
  while (std::getline(stream, line))
  {
    std::istringstream row(line);
    std::string first;
    std::string second;
    std::getline(row, first, ',');
    std::getline(row, second, ',');
    rows.emplace_back(std::stod(first), std::stod(second));
  }
  return rows;
}

void
expectPointwiseEqual(const std::vector<std::pair<Real, Real>> & expected,
                     const std::vector<Real> & values,
                     const std::vector<Real> & abscissa)
{
  ASSERT_EQ(expected.size(), values.size());
  ASSERT_EQ(expected.size(), abscissa.size());
  for (std::size_t i = 0; i < expected.size(); ++i)
  {
    EXPECT_DOUBLE_EQ(expected[i].first, values[i]);
    EXPECT_DOUBLE_EQ(expected[i].second, abscissa[i]);
  }
}

struct SyntheticFiles
{
  std::filesystem::path directory;
  std::string compression_hardening;
  std::string compression_damage;
  std::string tension_stiffening;
  std::string tension_damage;
};

void
writeFile(const std::filesystem::path & path, const std::string & content)
{
  std::ofstream stream(path);
  ASSERT_TRUE(stream.good());
  stream << content;
}

SyntheticFiles
syntheticFiles(const std::string & suffix)
{
  SyntheticFiles files;
  files.directory = std::filesystem::temp_directory_path() / ("cdp_material_table_" + suffix);
  std::filesystem::create_directories(files.directory);
  files.compression_hardening = (files.directory / "compression_hardening.csv").string();
  files.compression_damage = (files.directory / "compression_damage.csv").string();
  files.tension_stiffening = (files.directory / "tension_stiffening.csv").string();
  files.tension_damage = (files.directory / "tension_damage.csv").string();

  writeFile(files.compression_hardening, "stress_pa,inelastic_strain\n100,0\n80,0.001\n60,0.002\n");
  writeFile(files.compression_damage, "damage_c,inelastic_strain\n0,0\n0.1,0.001\n0.2,0.002\n");
  writeFile(files.tension_stiffening, "stress_pa,cracking_strain\n10,0\n8,0.001\n6,0.002\n");
  writeFile(files.tension_damage, "damage_t,cracking_strain\n0,0\n0.1,0.001\n0.2,0.002\n");
  return files;
}
}

TEST(CDPMaterialTable, ReadsReferencePythonBundlePointwise)
{
  const auto table = referenceTable();
  EXPECT_EQ(table.stressPointCount(CDPMaterialTable::Branch::COMPRESSION), 55);
  EXPECT_EQ(table.damagePointCount(CDPMaterialTable::Branch::COMPRESSION), 55);
  EXPECT_EQ(table.stressPointCount(CDPMaterialTable::Branch::TENSION), 56);
  EXPECT_EQ(table.damagePointCount(CDPMaterialTable::Branch::TENSION), 56);

  const auto compression_zero = table.response(CDPMaterialTable::Branch::COMPRESSION, 0.0);
  EXPECT_DOUBLE_EQ(compression_zero.stress.value, 8597100.0);
  EXPECT_DOUBLE_EQ(compression_zero.damage.value, 0.0);
  EXPECT_EQ(table.equivalentPlasticStrain(CDPMaterialTable::Branch::COMPRESSION).size(), 55);
  EXPECT_EQ(table.equivalentPlasticStrain(CDPMaterialTable::Branch::TENSION).size(), 56);

  expectPointwiseEqual(readPythonBundleTable(data_dir + "compression_hardening.csv"),
                       table.stressValues(CDPMaterialTable::Branch::COMPRESSION),
                       table.stressAbscissa(CDPMaterialTable::Branch::COMPRESSION));
  expectPointwiseEqual(readPythonBundleTable(data_dir + "compression_damage.csv"),
                       table.damageValues(CDPMaterialTable::Branch::COMPRESSION),
                       table.damageAbscissa(CDPMaterialTable::Branch::COMPRESSION));
  expectPointwiseEqual(readPythonBundleTable(data_dir + "tension_stiffening.csv"),
                       table.stressValues(CDPMaterialTable::Branch::TENSION),
                       table.stressAbscissa(CDPMaterialTable::Branch::TENSION));
  expectPointwiseEqual(readPythonBundleTable(data_dir + "tension_damage.csv"),
                       table.damageValues(CDPMaterialTable::Branch::TENSION),
                       table.damageAbscissa(CDPMaterialTable::Branch::TENSION));
}

TEST(CDPMaterialTable, InterpolatesWithExplicitOneSidedDerivativesAndConstantEndpoints)
{
  const auto files = syntheticFiles("interpolation");
  const CDPMaterialTable table(files.compression_hardening,
                               files.compression_damage,
                               files.tension_stiffening,
                               files.tension_damage,
                               1.0e9);

  const auto middle = table.response(CDPMaterialTable::Branch::COMPRESSION, 0.0005);
  EXPECT_DOUBLE_EQ(middle.stress.value, 90.0);
  EXPECT_DOUBLE_EQ(middle.damage.value, 0.05);
  EXPECT_DOUBLE_EQ(middle.stress.left_derivative, -20000.0);
  EXPECT_DOUBLE_EQ(middle.stress.right_derivative, -20000.0);

  const auto knot = table.response(CDPMaterialTable::Branch::COMPRESSION, 0.001);
  EXPECT_DOUBLE_EQ(knot.stress.value, 80.0);
  EXPECT_DOUBLE_EQ(knot.stress.left_derivative, -20000.0);
  EXPECT_DOUBLE_EQ(knot.stress.right_derivative, -20000.0);

  const auto before = table.response(CDPMaterialTable::Branch::COMPRESSION, -1.0);
  EXPECT_DOUBLE_EQ(before.stress.value, 100.0);
  EXPECT_DOUBLE_EQ(before.stress.left_derivative, 0.0);
  EXPECT_DOUBLE_EQ(before.stress.right_derivative, 0.0);

  const auto after = table.response(CDPMaterialTable::Branch::COMPRESSION, 1.0);
  EXPECT_DOUBLE_EQ(after.stress.value, 60.0);
  EXPECT_DOUBLE_EQ(after.stress.left_derivative, 0.0);
  EXPECT_DOUBLE_EQ(after.stress.right_derivative, 0.0);
}

TEST(CDPMaterialTable, RejectsNonIncreasingAbscissa)
{
  const auto files = syntheticFiles("non_increasing");
  writeFile(files.compression_hardening, "stress_pa,inelastic_strain\n100,0\n80,0\n");
  EXPECT_THROW(CDPMaterialTable(files.compression_hardening,
                                files.compression_damage,
                                files.tension_stiffening,
                                files.tension_damage,
                                1.0e9),
               std::invalid_argument);
}

TEST(CDPMaterialTable, RejectsDamageAtOrAboveOne)
{
  const auto files = syntheticFiles("damage_limit");
  writeFile(files.tension_damage, "damage_t,cracking_strain\n0,0\n1,0.001\n");
  EXPECT_THROW(CDPMaterialTable(files.compression_hardening,
                                files.compression_damage,
                                files.tension_stiffening,
                                files.tension_damage,
                                1.0e9),
               std::invalid_argument);
}

TEST(CDPMaterialTable, RejectsNegativeEquivalentPlasticStrain)
{
  const auto files = syntheticFiles("negative_plastic");
  writeFile(files.compression_hardening,
            "stress_pa,inelastic_strain\n100000000,0\n100000000,0.001\n");
  writeFile(files.compression_damage, "damage_c,inelastic_strain\n0,0\n0.9,0.001\n");
  EXPECT_THROW(CDPMaterialTable(files.compression_hardening,
                                files.compression_damage,
                                files.tension_stiffening,
                                files.tension_damage,
                                3.0e10),
               std::invalid_argument);
}
