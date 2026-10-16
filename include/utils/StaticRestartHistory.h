#pragma once
#include <stdexcept>

namespace StaticRestartHistory
{
// Assign each history directly: shifting once leaves an old static ramp in older.
template <typename Vector>
void initialize(const Vector & current, Vector & old, Vector * older = nullptr)
{
  if (&current == &old || (older && (older == &current || older == &old)))
    throw std::invalid_argument("Static restart histories must be distinct vectors");
  old = current;
  if (older)
    *older = current;
}
}
