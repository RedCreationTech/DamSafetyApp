#pragma once

/** Scalar Duvaut-Lions update, shared by the candidate and its seeded AD checks.
 * The caller validates times and damage bounds. Recovery may decrease this
 * scalar; irreversible branch histories remain separate. mu=0 takes priority.
 */
namespace CombinedDamageViscosity
{
template <typename Scalar>
Scalar update(const Scalar & old_damage, const Scalar & target, double h, double mu)
{
  if (mu == 0.0)
    return target;
  if (h == 0.0)
    return old_damage;
  return (mu * old_damage + h * target) / (mu + h);
}
}
