#include "SeagullApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
SeagullApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

SeagullApp::SeagullApp(InputParameters parameters) : MooseApp(parameters)
{
  SeagullApp::registerAll(_factory, _action_factory, _syntax);
}

SeagullApp::~SeagullApp() {}

void
SeagullApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"SeagullApp"});
  Registry::registerActionsTo(af, {"SeagullApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
SeagullApp::registerApps()
{
  registerApp(SeagullApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
SeagullApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  SeagullApp::registerAll(f, af, s);
}
extern "C" void
SeagullApp__registerApps()
{
  SeagullApp::registerApps();
}
