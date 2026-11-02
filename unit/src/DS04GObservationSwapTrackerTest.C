#include "CDPObservationSwapTracker.h"
#include "CDPObservationContextGuard.h"
#include "gtest/gtest.h"
#include <stdexcept>

TEST(DS04GObservationSwapTracker, PreservesSeparateCandidatesAndResetsPerElement)
{
  CDPObservationSwapTracker tracker(2);
  tracker.completed(9, 10, 11);
  tracker.completed(9, 10, 12);
  tracker.completed(8, 20, 21);
  auto a = tracker.take(9);
  ASSERT_EQ(a.calls.size(), 2u);
  EXPECT_EQ(a.calls[0].evaluation, 10u);
  EXPECT_EQ(a.calls[0].call, 11u);
  EXPECT_EQ(a.calls[1].call, 12u);
  EXPECT_FALSE(a.overflow);
  EXPECT_TRUE(tracker.take(9).calls.empty());
  auto b = tracker.take(8);
  ASSERT_EQ(b.calls.size(), 1u);
  EXPECT_EQ(b.calls[0].call, 21u);
}

TEST(DS04GObservationSwapTracker, OverflowIsExplicitAndDoesNotInventCandidates)
{
  CDPObservationSwapTracker tracker(1);
  tracker.completed(9, 10, 11);
  tracker.completed(9, 10, 12);
  const auto batch = tracker.take(9);
  ASSERT_EQ(batch.calls.size(), 1u);
  EXPECT_EQ(batch.calls[0].call, 11u);
  EXPECT_TRUE(batch.overflow);
  EXPECT_FALSE(tracker.take(9).overflow);
}

TEST(DS04GObservationSwapTracker, RecorderContextRestoredAfterFailure)
{
  std::uint64_t element = 7, call = 9;
  unsigned int qp = 2;
  try
  {
    CDPObservationContextGuard guard(element, qp, call);
    element = 99; qp = 0; call = 0;
    throw std::runtime_error("simulated observer failure");
  }
  catch (const std::runtime_error &) {}
  EXPECT_EQ(element, 7u);
  EXPECT_EQ(qp, 2u);
  EXPECT_EQ(call, 9u);
}
