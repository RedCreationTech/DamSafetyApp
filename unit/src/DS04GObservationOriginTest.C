#include "CDPObservation.h"
#include "CDPObservationOrigin.h"
#include "gtest/gtest.h"
#include <stdexcept>

TEST(DS04GObservationOrigin, LimitAndUnresolvedFramesAreExplicit)
{
  unsigned int emitted = 0;
  const auto capture = [](void ** frames, int capacity)
  {
    EXPECT_EQ(capacity, 3);
    for (int i = 0; i < capacity; ++i)
      frames[i] = reinterpret_cast<void *>(static_cast<std::uintptr_t>(i + 1));
    return capacity;
  };
  const auto resolve = [](void * address, CDPObservationOrigin::Frame & frame)
  {
    if (reinterpret_cast<std::uintptr_t>(address) == 2) return false;
    frame.module = "fake"; frame.offset = 123;
    return true; // a missing dynamic symbol is allowed but reported by the real sink
  };
  const auto status = CDPObservationOrigin::sample(2, capture, resolve,
      [&](unsigned int, const CDPObservationOrigin::Frame &, bool) { ++emitted; });
  EXPECT_EQ(status.frames, 2u);
  EXPECT_EQ(status.unresolved, 1u);
  EXPECT_TRUE(status.truncated);
  EXPECT_FALSE(status.failed);
  EXPECT_EQ(emitted, 2u);
}

TEST(DS04GObservationOrigin, InvalidEmptyAndThrownCapturesFailClosed)
{
  auto resolve = [](void *, CDPObservationOrigin::Frame &) { return true; };
  auto emit = [](unsigned int, const CDPObservationOrigin::Frame &, bool) {};
  auto empty = [](void **, int) { return 0; };
  EXPECT_TRUE(CDPObservationOrigin::sample(0, empty, resolve, emit).failed);
  EXPECT_TRUE(CDPObservationOrigin::sample(65, empty, resolve, emit).failed);
  EXPECT_TRUE(CDPObservationOrigin::sample(2, empty, resolve, emit).failed);
  auto throws = [](void **, int) -> int { throw std::runtime_error("test"); };
  EXPECT_TRUE(CDPObservationOrigin::sample(2, throws, resolve, emit).failed);
  auto one = [](void ** frames, int) { frames[0] = reinterpret_cast<void *>(1); return 1; };
  auto bad_resolver = [](void *, CDPObservationOrigin::Frame &) -> bool { throw std::runtime_error("resolve"); };
  EXPECT_TRUE(CDPObservationOrigin::sample(2, one, bad_resolver, emit).failed);
  auto bad_writer = [](unsigned int, const CDPObservationOrigin::Frame &, bool) { throw std::runtime_error("write"); };
  EXPECT_TRUE(CDPObservationOrigin::sample(2, one, resolve, bad_writer).failed);
}

TEST(DS04GObservationOrigin, DisabledAndOtherPointsDoNotInvokeSink)
{
  struct Sink : CDPObservation::Sink
  {
    unsigned int calls = 0;
    bool enterPoint(std::uint64_t e, unsigned int qp, bool) override { return e == 7 && qp == 0; }
    void traceOrigin() noexcept override { ++calls; }
    void write(const char *, const char *, const char *, const std::vector<double> &,
               const std::string &) noexcept override {}
  } sink;
  CDPObservation::configured = nullptr;
  CDPObservation::traceOrigin();
  EXPECT_EQ(sink.calls, 0u);
  CDPObservation::configured = &sink;
  { CDPObservation::PointScope wrong(8, 0); CDPObservation::traceOrigin(); }
  EXPECT_EQ(sink.calls, 0u);
  { CDPObservation::PointScope selected(7, 0); CDPObservation::traceOrigin(); }
  EXPECT_EQ(sink.calls, 1u);
  CDPObservation::configured = nullptr;
}
