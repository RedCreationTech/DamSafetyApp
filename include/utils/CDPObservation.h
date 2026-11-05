#pragma once

#include <array>
#include <cstdint>
#include <initializer_list>
#include <string>
#include <vector>

namespace CDPObservation
{
struct Context
{
  std::uint64_t local_solve = 0;
  unsigned int partition = 0;
  unsigned int substep = 0;
  unsigned int N = 0;
  unsigned int iteration = 0;
};

class Sink
{
public:
  virtual ~Sink() = default;
  virtual bool enterPoint(std::uint64_t element, unsigned int qp, bool new_call) noexcept = 0;
  virtual void materialCompleted() noexcept {}
  virtual void failed() noexcept {}
  virtual void write(const char * stream,
                     const char * event,
                     const char * name,
                     const std::vector<double> & values,
                     const std::string & note) noexcept = 0;
};

inline thread_local Sink * configured = nullptr;
inline thread_local Sink * current = nullptr;
inline thread_local Context context;
inline thread_local std::uint64_t sequence = 0;

inline bool
active()
{
  return current != nullptr;
}

inline void
materialCompleted() noexcept
{
  if (current)
    current->materialCompleted();
}

class PointScope
{
public:
  PointScope(std::uint64_t element, unsigned int qp, bool new_call = true) noexcept
    : _previous(current), _context(context)
  {
    current = configured && configured->enterPoint(element, qp, new_call) ? configured : nullptr;
    context = {};
  }

  ~PointScope()
  {
    current = _previous;
    context = _context;
  }

  PointScope(const PointScope &) = delete;
  PointScope & operator=(const PointScope &) = delete;

private:
  Sink * _previous;
  Context _context;
};

inline void
values(const char * stream,
       const char * event,
       const char * name,
       std::initializer_list<double> value,
       const std::string & note = {}) noexcept
{
  try
  {
    if (current)
      current->write(stream, event, name, std::vector<double>(value), note);
  }
  catch (...)
  {
    if (current)
      current->failed();
  }
}

template <std::size_t N>
void
array(const char * stream,
      const char * event,
      const char * name,
      const std::array<double, N> & value) noexcept
{
  try
  {
    if (current)
      current->write(stream, event, name, {value.begin(), value.end()}, {});
  }
  catch (...)
  {
    if (current)
      current->failed();
  }
}

template <std::size_t N, std::size_t M>
void
matrix(const char * stream,
       const char * event,
       const char * name,
       const std::array<std::array<double, M>, N> & value) noexcept
{
  try
  {
    if (!current)
      return;
    std::vector<double> flattened;
    flattened.reserve(N * M);
    for (const auto & row : value)
      flattened.insert(flattened.end(), row.begin(), row.end());
    current->write(stream, event, name, flattened, {});
  }
  catch (...)
  {
    if (current)
      current->failed();
  }
}
}
