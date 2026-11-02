#pragma once

#include <algorithm>

/** State transitions for a failure-triggered hold after an accepted recovery step.
 *
 * The caller owns and persists the two state variables. The hold counts only
 * successfully accepted steps AFTER the successful recovery step.
 */
namespace AqueductRetryCooldown
{
inline void failure(unsigned int requested, unsigned int & remaining, bool & recovering)
{
  remaining = requested;
  recovering = true;
}

inline void accepted(unsigned int & remaining, bool & recovering)
{
  if (recovering)
    recovering = false;
  else if (remaining)
    --remaining;
}

inline double limit(double proposed, double accepted_dt, unsigned int remaining)
{
  return remaining ? std::min(proposed, accepted_dt) : proposed;
}
}
