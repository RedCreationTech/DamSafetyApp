#include "CDPObservation.h"
#include "CDPObservationContextGuard.h"
#include "CDPObservationSwapTracker.h"

#include "gtest/gtest.h"

#include <array>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
struct ObservationSink : CDPObservation::Sink
{
  bool enterPoint(std::uint64_t element, unsigned int qp, bool) noexcept override
  {
    return element == 7 && qp == 0;
  }

  void write(const char *,
             const char * event,
             const char * field,
             const std::vector<double> & values,
             const std::string &) noexcept override
  {
    try
    {
      captured[std::string(event) + ":" + field] = values;
    }
    catch (...)
    {
    }
  }

  std::map<std::string, std::vector<double>> captured;
};
}

TEST(DS04GAcceptedOutputObservation, DisabledScopeDoesNotWrite)
{
  ObservationSink sink;
  CDPObservation::configured = nullptr;
  CDPObservation::current = nullptr;
  {
    CDPObservation::PointScope scope(7, 0);
    EXPECT_FALSE(CDPObservation::active());
    CDPObservation::values("qp-updates", "material_enter", "dt", {0.1});
  }
  EXPECT_TRUE(sink.captured.empty());
}

TEST(DS04GAcceptedOutputObservation, RotationMatrixKeepsAllNineRowMajorEntries)
{
  ObservationSink sink;
  CDPObservation::configured = &sink;
  {
    CDPObservation::PointScope scope(7, 0);
    const std::array<std::array<double, 3>, 3> rotation = {
        {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}}};
    CDPObservation::matrix("qp-updates", "material_enter", "rotation_row_major", rotation);
  }
  CDPObservation::configured = nullptr;
  const auto & values = sink.captured["material_enter:rotation_row_major"];
  ASSERT_EQ(values.size(), 9u);
  for (std::size_t i = 0; i < values.size(); ++i)
    EXPECT_EQ(values[i], static_cast<double>(i + 1));
}

TEST(DS04GAcceptedOutputObservation, PointAndFEContextsRestoreAfterExceptions)
{
  ObservationSink sink;
  CDPObservation::configured = &sink;
  CDPObservation::context.partition = 4;
  try
  {
    CDPObservation::PointScope scope(7, 0);
    CDPObservation::context.partition = 9;
    throw std::runtime_error("observer test");
  }
  catch (const std::runtime_error &)
  {
  }
  EXPECT_EQ(CDPObservation::context.partition, 4u);
  EXPECT_EQ(CDPObservation::current, nullptr);
  CDPObservation::configured = nullptr;
  CDPObservation::context = {};

  std::uint64_t evaluation = 0;
  const char * context = "execute";
  int iteration = -1;
  try
  {
    CDPObservationFEContextGuard guard(evaluation, context, iteration);
    evaluation = 12;
    context = "residual";
    iteration = 3;
    throw std::runtime_error("residual test");
  }
  catch (const std::runtime_error &)
  {
  }
  EXPECT_EQ(evaluation, 0u);
  EXPECT_STREQ(context, "execute");
  EXPECT_EQ(iteration, -1);
}

TEST(DS04GAcceptedOutputObservation, SwapCandidatesCarryAttemptExecutionAndCall)
{
  CDPObservationSwapTracker tracker(2);
  tracker.completed(7, 3, 5, 8, 13);
  tracker.completed(7, 4, 5, 9, 14);
  const auto batch = tracker.take(7);
  ASSERT_EQ(batch.calls.size(), 2u);
  EXPECT_EQ(batch.calls[0].attempt, 3u);
  EXPECT_EQ(batch.calls[0].execution, 5u);
  EXPECT_EQ(batch.calls[0].evaluation, 8u);
  EXPECT_EQ(batch.calls[0].call, 13u);
  EXPECT_EQ(batch.calls[1].call, 14u);
  EXPECT_FALSE(batch.overflow);
  EXPECT_TRUE(tracker.take(7).calls.empty());
}
