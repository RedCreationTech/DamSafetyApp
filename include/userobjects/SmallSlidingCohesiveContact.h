#pragma once
#include "MortarUserObject.h"
#include "AbaqusCohesiveLaw.h"
#include <unordered_map>

class SmallSlidingCohesiveContact : public MortarUserObject
{
public:
  static InputParameters validParams();
  SmallSlidingCohesiveContact(const InputParameters & parameters);
  void initialSetup() override;
  void initialize() override;
  void execute() override;
  void finalize() override;
  void reinit() override;
  const auto & weightedGaps() const { return _gaps; }
  const ADVariableValue & traction(unsigned int d) const { return _traction[d]; }
  Real nodeValue(const Node * node, unsigned int quantity) const;
protected:
  AbaqusCohesiveLaw::Result<ADReal> evaluateNode(const Node * node);
  std::array<Point,3> basis(const Node * node) const;
  const std::array<const ADVariableValue *,3> _secondary, _primary;
  const std::array<const MooseVariable *,3> _displacement_vars;
  const MooseVariable * _pressure_var;
  SystemBase & _system;
  const VariableTestValue * _shape=nullptr;
  const MooseArray<Real> & _coord;
  const AbaqusCohesiveLaw::Parameters _parameters;
  const Real _bond_tolerance;
  const unsigned int _expected_nodes;
  bool _coverage_reported=false;
  std::unordered_map<const DofObject *,std::pair<ADReal,Real>> _gaps;
  std::unordered_map<const DofObject *,ADRealVectorValue> _jumps, _reference_normals;
  std::unordered_map<const DofObject *,Real> _initial_gaps;
  std::unordered_map<dof_id_type,AbaqusCohesiveLaw::State> & _history;
  std::array<ADVariableValue,3> _traction;
};
