#include "IntegratedAcceleration.h"
#include <fstream>
#include <sstream>

registerMooseObject("DamSafetyApp", IntegratedAcceleration);
InputParameters IntegratedAcceleration::validParams()
{
  auto params = Function::validParams();
  params.addClassDescription("Analytic U/V/A from a headerless two-column acceleration CSV; independent of time-step partition.");
  params.addRequiredParam<FileName>("data_file", "Headerless ASCII time,acceleration file in SI units");
  params.addRequiredParam<MooseEnum>("component", MooseEnum("displacement velocity acceleration"), "Requested motion component");
  params.addParam<Real>("initial_displacement", 0, "Displacement at the first table time");
  params.addParam<Real>("initial_velocity", 0, "Velocity at the first table time");
  return params;
}
IntegratedAcceleration::IntegratedAcceleration(const InputParameters & parameters)
  : Function(parameters), _component(getParam<MooseEnum>("component"))
{
  const auto file = getParam<FileName>("data_file");
  std::ifstream input(file);
  if (!input) paramError("data_file", "Cannot open ", file);
  std::vector<double> times, accelerations;
  std::string line;
  unsigned int row = 0;
  while (std::getline(input, line))
  {
    ++row;
    if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream stream(line);
    double t, a;
    std::string extra;
    if (!(stream >> t >> a) || (stream >> extra))
      paramError("data_file", "Expected exactly two numeric columns at row ", row);
    times.push_back(t); accelerations.push_back(a);
  }
  try
  {
    _motion = std::make_unique<IntegratedGroundMotion>(times, accelerations,
        getParam<Real>("initial_displacement"), getParam<Real>("initial_velocity"));
  }
  catch (const std::exception & error) { paramError("data_file", error.what()); }
}
Real IntegratedAcceleration::value(Real t, const Point &) const
{
  const auto state = _motion->at(t);
  if (_component == "displacement") return state.displacement;
  if (_component == "velocity") return state.velocity;
  return state.acceleration;
}
Real IntegratedAcceleration::timeDerivative(Real t, const Point &) const
{
  const auto state = _motion->at(t);
  if (_component == "displacement") return state.velocity;
  if (_component == "velocity") return state.acceleration;
  return state.jerk;
}
