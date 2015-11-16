//
// $Id: B4dActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// ile B4dActionInitialization.cc
/// rief Implementation of the B4dActionInitialization class

#include "ExN01ActionInitialization.hh"
#include "ExN01PrimaryGeneratorAction.hh"
#include "ExN01RunAction.hh"
#include "ExN01EventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01ActionInitialization::ExN01ActionInitialization(UWUW *uwuw_workflow_data)
  : G4VUserActionInitialization()
{
  //  workflow_data = UWUW(filename);
  workflow_data = uwuw_workflow_data;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01ActionInitialization::~ExN01ActionInitialization()
{
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01ActionInitialization::BuildForMaster() const
{
  SetUserAction(new ExN01RunAction(workflow_data));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01ActionInitialization::Build() const
{
  SetUserAction(new ExN01PrimaryGeneratorAction);
  SetUserAction(new ExN01RunAction(workflow_data));
  SetUserAction(new ExN01EventAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
