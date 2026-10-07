#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

/** Stateless analytic integrals of a continuous, piecewise-linear acceleration.
 * No accepted/trial time-step history is stored. Queries outside the table fail
 * rather than silently extrapolating an earthquake record.
 */
class IntegratedGroundMotion
{
public:
  struct State { double displacement; double velocity; double acceleration; double jerk; };
  IntegratedGroundMotion(const std::vector<double> & times,
                         const std::vector<double> & accelerations,
                         double initial_displacement = 0, double initial_velocity = 0)
    : _times(times), _accelerations(accelerations),
      _displacements(times.size(), initial_displacement), _velocities(times.size(), initial_velocity)
  {
    if (times.size() < 2 || times.size() != accelerations.size() ||
        !std::isfinite(initial_displacement) || !std::isfinite(initial_velocity))
      throw std::invalid_argument("ground motion requires at least two finite (time, acceleration) pairs");
    for (std::size_t i = 0; i < times.size(); ++i)
    {
      if (!std::isfinite(times[i]) || !std::isfinite(accelerations[i]) ||
          (i && times[i] <= times[i-1]))
        throw std::invalid_argument("ground motion times must be finite and strictly increasing");
      if (i)
      {
        const double h = times[i] - times[i-1];
        const double slope = (accelerations[i] - accelerations[i-1]) / h;
        _velocities[i] = _velocities[i-1] + accelerations[i-1] * h + slope * h*h / 2;
        _displacements[i] = _displacements[i-1] + _velocities[i-1] * h +
                            accelerations[i-1] * h*h / 2 + slope * h*h*h / 6;
      }
    }
  }
  State at(double time) const
  {
    if (!std::isfinite(time) || time < _times.front() || time > _times.back())
      throw std::out_of_range("ground motion query outside acceleration time table");
    const auto upper = std::upper_bound(_times.begin(), _times.end(), time);
    const std::size_t i = upper == _times.end() ? _times.size()-2 : upper-_times.begin()-1;
    const double h = time - _times[i];
    const double slope = (_accelerations[i+1]-_accelerations[i])/(_times[i+1]-_times[i]);
    return {_displacements[i]+_velocities[i]*h+_accelerations[i]*h*h/2+slope*h*h*h/6,
            _velocities[i]+_accelerations[i]*h+slope*h*h/2, _accelerations[i]+slope*h, slope};
  }
private:
  std::vector<double> _times, _accelerations, _displacements, _velocities;
};
