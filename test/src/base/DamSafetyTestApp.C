//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "DamSafetyTestApp.h"
#include "DamSafetyApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
DamSafetyTestApp::validParams()
{
  InputParameters params = DamSafetyApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

DamSafetyTestApp::DamSafetyTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  DamSafetyTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

DamSafetyTestApp::~DamSafetyTestApp() {}

void
DamSafetyTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  DamSafetyApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"DamSafetyTestApp"});
    Registry::registerActionsTo(af, {"DamSafetyTestApp"});
  }
}

void
DamSafetyTestApp::registerApps()
{
  registerApp(DamSafetyApp);
  registerApp(DamSafetyTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
DamSafetyTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  DamSafetyTestApp::registerAll(f, af, s);
}
extern "C" void
DamSafetyTestApp__registerApps()
{
  DamSafetyTestApp::registerApps();
}
