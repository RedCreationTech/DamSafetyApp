#include "StaticRestartHistory.h"
#include "gtest/gtest.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"
#include <vector>
#include <cmath>

namespace
{
double hourglass(double q, double old, double older, double h)
{
  const double alpha = -.05, zeta = .00113;
  return q * (1 + alpha + (1 + alpha) * zeta / h) -
         old * (alpha + (1 + 2 * alpha) * zeta / h) + older * alpha * zeta / h;
}
}

TEST(StaticRestartHistory, LoadedDisplacementAndThicknessRemainAtRestForVariableIntervals)
{
  const std::vector<double> current{-.002315, .00037, -.000029};
  std::vector<double> old{-.0022, .00030, -.000023}, older{-.0021, .00028, -.000021};
  StaticRestartHistory::initialize(current, old, &older);
  EXPECT_EQ(old, current);
  EXPECT_EQ(older, current);
  for (double h : {.01, .005, .000234375})
    EXPECT_NEAR(hourglass(current[0], old[0], older[0], h), current[0], 1e-17);
}

TEST(StaticRestartHistory, SingleFrameworkShiftDoesNotEraseStaticRamp)
{
  const std::vector<double> current{1.0};
  std::vector<double> old{.99}, older{.98};
  older = old;
  old = current;
  EXPECT_NEAR(hourglass(1., old[0], older[0], .01), 1.0000565, 1e-15);
  StaticRestartHistory::initialize(current, old, &older);
  EXPECT_NEAR(hourglass(1., old[0], older[0], .01), 1., 1e-15);
}

TEST(StaticRestartHistory, RelativeAccelerationConventionIsPreserved)
{
  // vel_x, vel_y, relative accel_x free/base, relative accel_y.
  const std::vector<double> current{0., 0., .022, 0., 0.};
  std::vector<double> old(current.size(), 77.);
  StaticRestartHistory::initialize(current, old);
  EXPECT_EQ(old, current);
  EXPECT_DOUBLE_EQ(old[2] - .022, 0.); // absolute free-node initial acceleration
}

TEST(StaticRestartHistory, RejectedTrialDoesNotEnterAcceptedMotionHistory)
{
  std::vector<double> current{.002}, old{.0019}, older{.0018};
  StaticRestartHistory::initialize(current, old, &older);
  // Model the first normal advance, then a discarded trial and its restore.
  older = old;
  old = current;
  current[0] = .003;
  current = old;
  EXPECT_DOUBLE_EQ(old[0], .002);
  EXPECT_DOUBLE_EQ(older[0], .002);
  EXPECT_NEAR(hourglass(current[0], old[0], older[0], .005), .002, 1e-17);
}

TEST(StaticRestartHistory, AliasingIsRejectedBeforeMutation)
{
  std::vector<double> current{1.}, old{.9};
  EXPECT_THROW(StaticRestartHistory::initialize(current, current, &old), std::invalid_argument);
  EXPECT_THROW(StaticRestartHistory::initialize(current, old, &old), std::invalid_argument);
  EXPECT_DOUBLE_EQ(current[0], 1.);
  EXPECT_DOUBLE_EQ(old[0], .9);
}

TEST(StaticRestartHistory, ActualLibMeshParallelToGhostedHistoryPreservesCurrent)
{
  libMesh::Parallel::Communicator comm(MPI_COMM_WORLD);
  const libMesh::numeric_index_type local = 3, size = local * comm.size();
  const auto first = local * comm.rank();
  std::vector<libMesh::numeric_index_type> ghosts;
  if (comm.size() > 1)
    ghosts.push_back(local * ((comm.rank() + 1) % comm.size()));
  auto current = libMesh::NumericVector<double>::build(comm);
  auto old = libMesh::NumericVector<double>::build(comm);
  auto older = libMesh::NumericVector<double>::build(comm);
  current->init(size, local, false, libMesh::PARALLEL);
  old->init(size, local, ghosts, false, libMesh::GHOSTED);
  older->init(size, local, ghosts, false, libMesh::GHOSTED);
  for (auto i = first; i < first + local; ++i)
  {
    current->set(i, .002315 * (i + 1));
    old->set(i, -.5);
    older->set(i, -.7);
  }
  current->close(); old->close(); older->close();
  StaticRestartHistory::initialize(*current, *old, older.get());
  old->close(); older->close();
  for (auto i = first; i < first + local; ++i)
  {
    EXPECT_DOUBLE_EQ((*current)(i), .002315 * (i + 1));
    EXPECT_DOUBLE_EQ((*old)(i), (*current)(i));
    EXPECT_DOUBLE_EQ((*older)(i), (*current)(i));
  }
  for (const auto i : ghosts)
  {
    EXPECT_DOUBLE_EQ((*old)(i), .002315 * (i + 1));
    EXPECT_DOUBLE_EQ((*older)(i), .002315 * (i + 1));
  }
}
