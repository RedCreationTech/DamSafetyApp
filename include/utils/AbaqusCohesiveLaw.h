#pragma once
#include "MooseTypes.h"
#include "ADReal.h"
#include <array>
#include <algorithm>
#include <cmath>

// The state is committed only after an accepted time step. Trial evaluations
// are pure functions of that committed state, including damage initiation.
namespace AbaqusCohesiveLaw
{
using State = std::array<Real, 9>;
struct Parameters
{
  std::array<Real,3> stiffness, strength;
  Real failure_increment, viscosity, friction, critical_slip;
};
template <typename T> struct Result
{
  std::array<T,3> traction;
  T damage;
  State state;
};
template <typename T>
Result<T> evaluate(const Parameters & p, const State & old,
                   const std::array<T,3> & jump, const T & pressure, Real dt, bool bonded=true)
{
  using std::sqrt; using std::abs;
  // Select the tensile generalized derivative at zero separation. This avoids
  // an artificial zero tangent for a closed, initially undamaged bond.
  auto positive=[](const T & x)->T { return x >= 0 ? x : T(0); };
  const T opening=positive(jump[0]);
  const T effective=sqrt(opening*opening+jump[1]*jump[1]+jump[2]*jump[2]);
  T onset=old[0];
  if (bonded && old[0]==0)
  {
    const T ratio=std::max(positive(p.stiffness[0]*jump[0])/p.strength[0],
                    std::max(abs(p.stiffness[1]*jump[1])/p.strength[1],
                             abs(p.stiffness[2]*jump[2])/p.strength[2]));
    if (ratio>=1)
    {
      // Interpolate the earliest MAXS crossing along this accepted-step trial
      // path; radial rescaling alone is wrong for a nonproportional increment.
      T alpha=1;
      for (unsigned int d=0;d<3;++d)
        for (int sign : {-1,1})
        {
          if (d==0 && sign==-1) continue;
          const Real threshold=sign*p.strength[d]/p.stiffness[d];
          const T delta=jump[d]-old[2+d];
          if ((sign*jump[d]>=sign*threshold) && (sign*old[2+d]<sign*threshold))
            alpha=std::min(alpha,(threshold-old[2+d])/delta);
        }
      std::array<T,3> first;
      for (unsigned int d=0;d<3;++d) first[d]=old[2+d]+alpha*(jump[d]-old[2+d]);
      first[0]=positive(first[0]);
      onset=sqrt(first[0]*first[0]+first[1]*first[1]+first[2]*first[2]);
    }
  }
  const T maximum=std::max(T(old[1]),effective);
  T damage_target=old[8];
  if (!bonded) damage_target=1;
  else if (onset>0 && maximum>onset)
  {
    const T final=onset+p.failure_increment;
    const T candidate=maximum>=final ? T(1) :
        final*(maximum-onset)/(maximum*p.failure_increment);
    damage_target=std::max(damage_target,candidate);
  }
  const T damage=bonded ? (dt*damage_target+p.viscosity*old[7])/(dt+p.viscosity) : T(1);
  Result<T> out;
  out.damage=damage;
  out.traction={-(1-damage)*p.stiffness[0]*opening,
                -(1-damage)*p.stiffness[1]*jump[1],
                -(1-damage)*p.stiffness[2]*jump[2]};
  std::array<T,2> elastic{0,0};
  const T contact=positive(pressure);
  out.traction[0]+=pressure;
  if (contact>0 && p.friction>0)
  {
    for (unsigned int d=0;d<2;++d) elastic[d]=old[5+d]+jump[1+d]-old[3+d];
    const T norm=sqrt(elastic[0]*elastic[0]+elastic[1]*elastic[1]);
    if (norm>p.critical_slip)
      for (auto & value : elastic) value*=p.critical_slip/norm;
    for (unsigned int d=0;d<2;++d)
      out.traction[d+1]-=damage*p.friction*contact*elastic[d]/p.critical_slip;
  }
  out.state=old;
  out.state[0]=MetaPhysicL::raw_value(onset); out.state[1]=MetaPhysicL::raw_value(maximum);
  for(unsigned int d=0;d<3;++d) out.state[2+d]=MetaPhysicL::raw_value(jump[d]);
  for(unsigned int d=0;d<2;++d) out.state[5+d]=MetaPhysicL::raw_value(elastic[d]);
  out.state[7]=MetaPhysicL::raw_value(damage); out.state[8]=MetaPhysicL::raw_value(damage_target);
  return out;
}
}
