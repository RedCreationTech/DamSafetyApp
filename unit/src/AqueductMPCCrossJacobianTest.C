#include "AqueductMPCCrossJacobian.h"
#include "AqueductMPCAssemblyStorage.h"
#include "gtest/gtest.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/parallel.h"
#include <array>
#include <tuple>

class AqueductMPCCrossJacobianTest : public ::testing::Test
{
protected:
  libMesh::Parallel::Communicator comm{MPI_COMM_WORLD};
  libMesh::ReplicatedMesh mesh{comm, 1};
  libMesh::EquationSystems equations{mesh};
  libMesh::ExplicitSystem * system = nullptr;

  void SetUp() override
  {
    libMesh::MeshTools::Generation::build_line(mesh, 9, 0., 1., libMesh::EDGE2);
    system = &equations.add_system<libMesh::ExplicitSystem>("cross_block_contract");
    system->add_variable("u", libMesh::FIRST);
    equations.init();
  }
};

TEST_F(AqueductMPCCrossJacobianTest, MatchesIndependentResidualDerivativeWithDuplicateDofs)
{
  const std::vector<double> c{1., -.3, -.7, .1};
  const std::array<double, 4> scales{2., .2, 1., 3.};
  const double lm_scale = .7, compliance = 1. / 605000.;
  const std::vector<libMesh::dof_id_type> ids{2, 3, 2, 4, 5};
  AqueductMPCAssemblyStorage storage(c);
  const auto & dense = storage.jacobian(scales, lm_scale, compliance);
  AqueductMPCCrossJacobian cross(c.size());
  std::array<std::array<double, 10>, 10> actual{};
  cross.assemble(dense, ids, system->get_dof_map(),
                 [&](const auto & block, const auto & rows, const auto & cols) {
                   for (std::size_t i = 0; i < rows.size(); ++i)
                     for (std::size_t j = 0; j < cols.size(); ++j)
                       actual[rows[i]][cols[j]] += block(i, j);
                 });
  auto residual = [&](const std::array<double, 10> & q) {
    std::array<double, 10> r{};
    double gap = -compliance * q[ids.back()];
    for (std::size_t i = 0; i < c.size(); ++i)
    {
      r[ids[i]] += c[i] * scales[i] * q[ids.back()];
      gap += c[i] * q[ids[i]];
    }
    r[ids.back()] += gap * lm_scale;
    return r;
  };
  std::array<double, 10> q{};
  q[2] = .01; q[3] = -.02; q[4] = .004; q[5] = 4.;
  const double h = 1e-6;
  for (std::size_t j = 0; j < q.size(); ++j)
  {
    auto plus = q, minus = q;
    plus[j] += h; minus[j] -= h;
    const auto rp = residual(plus), rm = residual(minus);
    for (std::size_t i = 0; i < q.size(); ++i)
      EXPECT_NEAR(actual[i][j], (rp[i] - rm[i]) / (2 * h), 1e-9);
  }
}

TEST_F(AqueductMPCCrossJacobianTest, PreservesNonzeroInsertionOrderAndTinyCoefficients)
{
  const std::vector<libMesh::dof_id_type> ids{2, 3, 4, 5};
  AqueductMPCAssemblyStorage storage({1., 1e-200, -.7});
  const std::array<double, 3> scales{2., .2, 3.};
  const auto & dense = storage.jacobian(scales, .7, 1. / 21965000.);
  const auto original = dense;
  using Entry = std::tuple<libMesh::dof_id_type, libMesh::dof_id_type, double>;
  std::vector<Entry> expected, actual;
  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 4; ++j)
      if (dense(i, j) != 0.) expected.emplace_back(ids[i], ids[j], dense(i, j));
  AqueductMPCCrossJacobian cross(3);
  std::size_t calls = 0, slots = 0;
  cross.assemble(dense, ids, system->get_dof_map(),
                 [&](const auto & block, const auto & rows, const auto & cols) {
                   ++calls; slots += block.m() * block.n();
                   for (std::size_t i = 0; i < rows.size(); ++i)
                     for (std::size_t j = 0; j < cols.size(); ++j)
                       if (block(i, j) != 0.)
                         actual.emplace_back(rows[i], cols[j], block(i, j));
                 });
  EXPECT_EQ(actual, expected);
  EXPECT_EQ(calls, 3u);
  EXPECT_EQ(slots, 7u);
  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 4; ++j)
      EXPECT_DOUBLE_EQ(dense(i, j), original(i, j));
}

TEST_F(AqueductMPCCrossJacobianTest, RefreshesDofsScalesAndComplianceOnRetry)
{
  AqueductMPCAssemblyStorage storage({1., -.2, -.8});
  AqueductMPCCrossJacobian cross(3);
  std::array<double, 3> scales{1., 2., 3.};
  for (unsigned round = 0; round < 3; ++round)
  {
    const std::vector<libMesh::dof_id_type> ids =
        round == 1 ? std::vector<libMesh::dof_id_type>{6, 5, 4, 3}
                   : std::vector<libMesh::dof_id_type>{2, 3, 4, 5};
    scales[1] = round == 1 ? 7. : 2.;
    const double compliance = round == 1 ? .25 : 0.;
    const auto & dense = storage.jacobian(scales, 4., compliance);
    std::size_t calls = 0;
    cross.assemble(dense, ids, system->get_dof_map(),
                   [&](const auto & block, const auto & rows, const auto & cols) {
                     if (calls == 0)
                     {
                       EXPECT_EQ(rows, (std::vector<libMesh::dof_id_type>(ids.begin(), ids.end()-1)));
                       EXPECT_EQ(cols, (std::vector<libMesh::dof_id_type>{ids.back()}));
                       EXPECT_DOUBLE_EQ(block(1, 0), -.2 * scales[1]);
                     }
                     if (calls == 2) EXPECT_DOUBLE_EQ(block(0, 0), -1.);
                     ++calls;
                   });
    EXPECT_EQ(calls, compliance ? 3u : 2u);
  }
}

TEST_F(AqueductMPCCrossJacobianTest, RealLibMeshConstrainedDofUsesUnchangedDensePath)
{
  libMesh::DofConstraintRow relation;
  relation[6] = .25; relation[7] = .75;
  system->get_dof_map().add_constraint_row(3, relation);
  const std::vector<libMesh::dof_id_type> ids{2, 3, 4, 5};
  AqueductMPCAssemblyStorage storage({1., -.2, -.8});
  const auto & dense = storage.jacobian(std::array<double, 3>{1., 2., 3.}, 4., .25);
  AqueductMPCCrossJacobian cross(3);
  std::size_t calls = 0;
  cross.assemble(dense, ids, system->get_dof_map(),
                 [&](const auto & block, const auto & rows, const auto & cols) {
                   ++calls;
                   EXPECT_EQ(&block, &dense);
                   EXPECT_EQ(rows, ids); EXPECT_EQ(cols, ids);
                   auto expected = dense, actual = block;
                   auto er = ids, ec = ids, ar = rows, ac = cols;
                   system->get_dof_map().constrain_element_matrix(expected, er, ec, false);
                   system->get_dof_map().constrain_element_matrix(actual, ar, ac, false);
                   EXPECT_EQ(er, ar); EXPECT_EQ(ec, ac);
                   for (unsigned i = 0; i < expected.m(); ++i)
                     for (unsigned j = 0; j < expected.n(); ++j)
                       EXPECT_DOUBLE_EQ(actual(i, j), expected(i, j));
                 });
  EXPECT_EQ(calls, 1u);
}

TEST_F(AqueductMPCCrossJacobianTest, RejectsWrongDimensions)
{
  EXPECT_THROW(AqueductMPCCrossJacobian(0), std::invalid_argument);
  AqueductMPCCrossJacobian cross(3);
  libMesh::DenseMatrix<libMesh::Number> matrix(4, 4);
  auto unused = [](const auto &, const auto &, const auto &) {};
  EXPECT_THROW(cross.assemble(matrix, {2, 3, 4}, system->get_dof_map(), unused),
               std::invalid_argument);
  matrix.resize(3, 4);
  EXPECT_THROW(cross.assemble(matrix, {2, 3, 4, 5}, system->get_dof_map(), unused),
               std::invalid_argument);
}
