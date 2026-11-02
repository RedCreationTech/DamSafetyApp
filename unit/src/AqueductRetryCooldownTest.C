#include "AqueductRetryCooldown.h"
#include "gtest/gtest.h"

TEST(AqueductRetryCooldown, RecoveryAndThreeAdditionalAcceptedStepsThenGrowth)
{
  unsigned int remaining = 0;
  bool recovering = false;
  AqueductRetryCooldown::failure(3, remaining, recovering);
  ASSERT_TRUE(recovering);
  ASSERT_EQ(remaining, 3u);

  // Recovery acceptance is not one of the three subsequent hold steps.
  AqueductRetryCooldown::accepted(remaining, recovering);
  EXPECT_FALSE(recovering);
  for (unsigned int i = 0; i < 3; ++i)
  {
    EXPECT_DOUBLE_EQ(AqueductRetryCooldown::limit(.0025, .00125, remaining), .00125);
    AqueductRetryCooldown::accepted(remaining, recovering);
  }
  EXPECT_EQ(remaining, 0u);
  EXPECT_DOUBLE_EQ(AqueductRetryCooldown::limit(.0025, .00125, remaining), .0025);
}

TEST(AqueductRetryCooldown, RepeatedFailureRestartsHoldWithoutUnsignedOverflow)
{
  unsigned int remaining = 0;
  bool recovering = false;
  AqueductRetryCooldown::failure(2, remaining, recovering);
  AqueductRetryCooldown::accepted(remaining, recovering);
  AqueductRetryCooldown::accepted(remaining, recovering);
  ASSERT_EQ(remaining, 1u);
  AqueductRetryCooldown::failure(2, remaining, recovering);
  ASSERT_EQ(remaining, 2u);
  ASSERT_TRUE(recovering);
  AqueductRetryCooldown::accepted(remaining, recovering);
  ASSERT_EQ(remaining, 2u);
  AqueductRetryCooldown::accepted(remaining, recovering);
  AqueductRetryCooldown::accepted(remaining, recovering);
  AqueductRetryCooldown::accepted(remaining, recovering);
  EXPECT_EQ(remaining, 0u);
}

TEST(AqueductRetryCooldown, AShorterKnotOrEndpointIsNeverRaised)
{
  EXPECT_DOUBLE_EQ(AqueductRetryCooldown::limit(.0001, .00125, 3), .0001);
  // The stored recovery ceiling remains .00125 after a knot step of .0001.
  EXPECT_DOUBLE_EQ(AqueductRetryCooldown::limit(.0025, .00125, 2), .00125);
  EXPECT_DOUBLE_EQ(AqueductRetryCooldown::limit(.0001, .00125, 0), .0001);
  EXPECT_DOUBLE_EQ(AqueductRetryCooldown::limit(.0025, .00125, 0), .0025);
}

TEST(AqueductRetryCooldown, PersistedStateContinuesCountingAfterRecovery)
{
  unsigned int original_remaining = 0;
  bool original_recovering = false;
  AqueductRetryCooldown::failure(3, original_remaining, original_recovering);
  AqueductRetryCooldown::accepted(original_remaining, original_recovering);
  AqueductRetryCooldown::accepted(original_remaining, original_recovering);
  // Model restoring both values from a checkpoint; actual MOOSE recovery is a separate test.
  unsigned int restored_remaining = original_remaining;
  bool restored_recovering = original_recovering;
  EXPECT_EQ(restored_remaining, 2u);
  EXPECT_DOUBLE_EQ(AqueductRetryCooldown::limit(.0025, .00125, restored_remaining), .00125);
  AqueductRetryCooldown::accepted(restored_remaining, restored_recovering);
  AqueductRetryCooldown::accepted(restored_remaining, restored_recovering);
  EXPECT_DOUBLE_EQ(AqueductRetryCooldown::limit(.0025, .00125, restored_remaining), .0025);
}
