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

  // set mandatory initialization classes
  //
  //G4VUserDetectorConstruction* detector = new ExN01DetectorConstruction(uwuw_file);
  runManager->SetUserInitialization(new ExN01DetectorConstruction(uwuw_file));
  //  runManager->SetUserInitialization(detector);

  G4PhysListFactory *physListFactory = new G4PhysListFactory();
  G4VUserPhysicsList *physicsList =
    physListFactory->GetReferencePhysList("QGSP_BIC_HP");
  runManager->SetUserInitialization(physicsList);


  // set mandatory user action class
  //
  ExN01ActionInitialization* actionInitialization = new ExN01ActionInitialization(uwuw_file);
  runManager->SetUserInitialization(actionInitialization);

  //  G4VUserPrimaryGeneratorAction* gen_action = new ExN01PrimaryGeneratorAction;
  //runManager->SetUserAction(gen_action);


  //  G4VUserSteppingAction* step_action = new ExN01SteppingAction;
  // runManager->SetUserAction(step_action);

  // Initialize G4 kernel
  //
  runManager->Initialize();


  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();


  // Get the pointer to the UI manager and set verbosities
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4UIExecutive* UI = new G4UIExecutive(argc, argv);
  UImanager->ApplyCommand("/control/execute vis.mac");
  UI->SessionStart();

  delete UI;

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

  return 0;


  UImanager->ApplyCommand("/run/verbose 1");
  UImanager->ApplyCommand("/event/verbose 1");
  UImanager->ApplyCommand("/tracking/verbose 1");



  // Start a run
  //
  G4int numberOfEvent = 1e6;
  runManager->BeamOn(numberOfEvent);

  // Job termination
  //
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  //
  delete runManager;

  return 0;
}
