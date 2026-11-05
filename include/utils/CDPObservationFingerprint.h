#pragma once

#include <cstdint>
#include <cstring>

class CDPObservationFingerprint
{
public:
  CDPObservationFingerprint(std::uint64_t first, std::uint64_t last)
  {
    addWord(first);
    addWord(last);
  }

  void add(double value)
  {
    static_assert(sizeof(double) == sizeof(std::uint64_t), "Requires 64-bit double");
    std::uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    addWord(bits);
  }

  std::uint64_t value() const { return _value; }

private:
  void addWord(std::uint64_t word)
  {
    for (unsigned int i = 0; i < 8; ++i)
    {
      _value ^= (word >> (8 * i)) & 0xff;
      _value *= 1099511628211ULL;
    }
  }

  std::uint64_t _value = 14695981039346656037ULL;
};
