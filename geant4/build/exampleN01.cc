//
//
// --------------------------------------------------------------
//      GEANT 4 - exampleN01
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4Timer.hh"

#include "ExN01DetectorConstruction.hh"
#include "ExN01PhysicsList.hh"
#include "ExN01PrimaryGeneratorAction.hh"
#include "ExN01ActionInitialization.hh"

#include "ExN01UserScoreWriter.hh"

#include "G4PhysListFactory.hh"
#include "G4ScoringManager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#ifndef uwuw_hpp
#define uwuw_hpp 1
#include "uwuw.hpp"
#endif

int main(int argc, char* argv[])
{
  G4Timer Timer;
  Timer.Start();

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // Activate command-based scorer
  G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
  scManager->SetVerboseLevel(1);
  scManager->SetScoreWriter(new ExN01UserScoreWriter());

  std::string uwuw_file(argv[1]); // file containing data & uwuw

  // Activate UI-command base scorer
  // load the UWUW data
  UWUW *workflow_data = new UWUW(uwuw_file);

  // setup detectors and scores
  runManager->SetUserInitialization(new ExN01DetectorConstruction(workflow_data));

  G4PhysListFactory *physListFactory = new G4PhysListFactory();
  G4VUserPhysicsList *physicsList =
      physListFactory->GetReferencePhysList("QGSP_BIC_HP");
  runManager->SetUserInitialization(physicsList);


  // set mandatory user action class
  //
  ExN01ActionInitialization* actionInitialization = new ExN01ActionInitialization(workflow_data);
  runManager->SetUserInitialization(actionInitialization);

  //  G4VUserPrimaryGeneratorAction* gen_action = new ExN01PrimaryGeneratorAction;
  //runManager->SetUserAction(gen_action);


  //  G4VUserSteppingAction* step_action = new ExN01SteppingAction;
  // runManager->SetUserAction(step_action);

  // Initialize G4 kernel
  //
  runManager->Initialize();

  // Get the pointer to the UI manager and set verbosities
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // batch mode
  if( argc > 2 ) {
    G4String command = "/control/execute ";
    std::string filename(argv[2]);
    G4UIExecutive* ui = new G4UIExecutive(argc,argv, "tcsh" );
    UImanager->ApplyCommand(command+filename);
    ui->SessionStart();
    delete ui;
  } else {
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    G4UIExecutive* UI = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute vis.mac");

    UI->SessionStart();
    delete visManager;
    delete UI;
  }

  // stop the timer
  Timer.Stop();
  G4cout << G4endl;
  G4cout << "******************************************";
  G4cout << G4endl;
  G4cout << "Total Real Elapsed Time is: "<< Timer.GetRealElapsed();
  G4cout << G4endl;
  G4cout << "Total System Elapsed Time: " << Timer.GetSystemElapsed();
  G4cout << G4endl;
  G4cout << "Total GetUserElapsed Time: " << Timer.GetUserElapsed();
  G4cout << G4endl;
  G4cout << "******************************************";
  G4cout << G4endl;

  delete runManager;

  return 0;
}
