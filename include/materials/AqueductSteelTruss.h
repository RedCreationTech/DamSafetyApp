#pragma once
#include "TrussMaterial.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace AqueductTruss
{
struct State { double stress,plastic,equivalent,tangent; };
inline State update(double strain,double plastic_old,double equivalent_old,double E,double fy,double fu,double limit)
{
  if(!(E>0 && fy>0 && fu>=fy && limit>0 && equivalent_old>=0) ||
     !std::isfinite(strain) || !std::isfinite(plastic_old) || !std::isfinite(equivalent_old) ||
     !std::isfinite(E) || !std::isfinite(fy) || !std::isfinite(fu) || !std::isfinite(limit))
    throw std::invalid_argument("invalid aqueduct truss material/state");
  const double H=(fu-fy)/limit,strength=fy+H*std::min(equivalent_old,limit);
  const double trial=E*(strain-plastic_old);
  if(std::abs(trial)<=strength)return {trial,plastic_old,equivalent_old,E};
  double slope=equivalent_old<limit?H:0;
  double increment=(std::abs(trial)-strength)/(E+slope);
  if(equivalent_old+increment>limit){increment=(std::abs(trial)-fu)/E;slope=0;}
  double plastic=plastic_old+std::copysign(increment,trial);
  return {trial-E*std::copysign(increment,trial),plastic,equivalent_old+increment,E*slope/(E+slope)};
}
}

// Source NLGEOM=NO: initial-axis strain/orientation, isotropic accumulated plasticity.
class AqueductSteelTruss : public TrussMaterial
{
public:
  static InputParameters validParams();
  AqueductSteelTruss(const InputParameters &);
protected:
  void initQpStatefulProperties() override;
  void computeQpStrain() override;
  void computeQpStress() override;
  Real _fy,_fu,_limit;
  MaterialProperty<Real> & _plastic;
  const MaterialProperty<Real> & _plastic_old;
  MaterialProperty<Real> & _equivalent;
  const MaterialProperty<Real> & _equivalent_old;
};
