#include "DamSafetyApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
DamSafetyApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

DamSafetyApp::DamSafetyApp(const InputParameters & parameters) : MooseApp(parameters)
{
  DamSafetyApp::registerAll(_factory, _action_factory, _syntax);
}

DamSafetyApp::~DamSafetyApp() {}

void
DamSafetyApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<DamSafetyApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"DamSafetyApp"});
  Registry::registerActionsTo(af, {"DamSafetyApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
DamSafetyApp::registerApps()
{
  registerApp(DamSafetyApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
DamSafetyApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  DamSafetyApp::registerAll(f, af, s);
}
extern "C" void
DamSafetyApp__registerApps()
{
  DamSafetyApp::registerApps();
}
