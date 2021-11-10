//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "SeagullTestApp.h"
#include "SeagullApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
SeagullTestApp::validParams()
{
  InputParameters params = SeagullApp::validParams();
  return params;
}

SeagullTestApp::SeagullTestApp(InputParameters parameters) : MooseApp(parameters)
{
  SeagullTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

SeagullTestApp::~SeagullTestApp() {}

void
SeagullTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  SeagullApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"SeagullTestApp"});
    Registry::registerActionsTo(af, {"SeagullTestApp"});
  }
}

void
SeagullTestApp::registerApps()
{
  registerApp(SeagullApp);
  registerApp(SeagullTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
SeagullTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  SeagullTestApp::registerAll(f, af, s);
}
extern "C" void
SeagullTestApp__registerApps()
{
  SeagullTestApp::registerApps();
}
