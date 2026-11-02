#pragma once
#include <array>
#include <cstdint>

/** Bounded, injectable stack collection; this is diagnostic evidence, not identity proof. */
namespace CDPObservationOrigin
{
constexpr unsigned int max_frames = 64;
struct Frame
{
  const char * module = nullptr;
  const char * symbol = nullptr;
  std::uintptr_t offset = 0;
};
struct Status
{
  unsigned int frames = 0;
  unsigned int unresolved = 0;
  bool truncated = false;
  bool failed = false;
};

template <typename Capture, typename Resolve, typename Emit>
Status sample(unsigned int limit, Capture capture, Resolve resolve, Emit emit) noexcept
{
  Status out;
  if (limit == 0 || limit > max_frames) { out.failed = true; return out; }
  try
  {
    // One extra slot detects truncation; never request more than this bound.
    std::array<void *, max_frames + 1> addresses{};
    const int n = capture(addresses.data(), static_cast<int>(limit + 1));
    if (n <= 0 || n > static_cast<int>(limit + 1)) { out.failed = true; return out; }
    out.truncated = n > static_cast<int>(limit);
    const auto used = out.truncated ? limit : static_cast<unsigned int>(n);
    for (unsigned int i = 0; i < used; ++i)
    {
      Frame frame;
      const bool resolved = resolve(addresses[i], frame);
      if (!resolved) ++out.unresolved;
      emit(i, frame, resolved);
      ++out.frames;
    }
  }
  catch (...) { out.failed = true; }
  return out;
}
} // namespace CDPObservationOrigin
