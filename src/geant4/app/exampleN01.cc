//
//
// --------------------------------------------------------------
//      GEANT 4 - exampleN01
// --------------------------------------------------------------

#include <filesystem>

#include "ExN01ActionInitialization.hh"
#include "ExN01DetectorConstruction.hh"
#include "ExN01PhysicsList.hh"
#include "ExN01PrimaryGeneratorAction.hh"
#include "ExN01UserScoreWriter.hh"
#include "G4PhysListFactory.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4ScoringManager.hh"
#include "G4Timer.hh"
#include "G4UImanager.hh"

#ifdef GEANT4_GT_10_6
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#else
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
#endif

#ifndef uwuw_hpp
#define uwuw_hpp 1
#include "uwuw.hpp"
#endif

#include "moab/ProgOptions.hpp" // prog options

int main(int argc, char* argv[]) {
  G4Timer Timer;
  Timer.Start();

  // Activate UI-command base scorer
  // load the UWUW data

  std::string dag_file = "";
  std::string bfield_file = "";
  std::string macro_file = "";
  
  ProgOptions po("DagGeant4: a DAG tool for Geant4");
  po.addOpt<std::string>("dagmc,d", "Path to DAGMC file to proccess", &dag_file);
  po.addOpt<std::string>("bfield,b", "Path to Bfield file to proccess", &bfield_file);
  po.addOpt<std::string>("macro,m", "Path to DAGMC file to proccess", &macro_file);

  po.addOptionHelpHeading("Options for loading files");

  // get the options
  po.parseCommandLine(argc, argv);
  UWUW *workflow_data = NULL;
  // check too see the dagmc file exists
  if(std::filesystem::exists(dag_file)){
    workflow_data = new UWUW(dag_file);
  } else {
    std::cerr << "Error: dagmc file does not exist" << std::endl;
    return -1;
  }

  // if the bfield string set - see if the file exists
  if(bfield_file.length() && !std::filesystem::exists(bfield_file)){
    std::cerr << "Error: bfield file does not exist" << std::endl;
    return -1;
  }

  // detect if we are in batch mode or not
  G4UIExecutive* ui = nullptr;
  if (macro_file.length() == 0) {
    ui = new G4UIExecutive(argc, argv);
  }

  //G4RunManager* runManager =
  //    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
  G4RunManager* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Tasking);

  // Activate command-based scorer
  G4ScoringManager* scManager = G4ScoringManager::GetScoringManager();
  scManager->SetVerboseLevel(20);
  scManager->SetScoreWriter(new ExN01UserScoreWriter());
 
  // setup detectors and scores
  runManager->SetUserInitialization(
      new ExN01DetectorConstruction(workflow_data, bfield_file));

  G4PhysListFactory* physListFactory = new G4PhysListFactory();
  G4VUserPhysicsList* physicsList =
      physListFactory->GetReferencePhysList("QGSP_BIC_AllHPT");
  runManager->SetUserInitialization(physicsList);

  // set mandatory user action class
  ExN01ActionInitialization* actionInitialization =
      new ExN01ActionInitialization(workflow_data);
  runManager->SetUserInitialization(actionInitialization);

  // lets get started
  //runManager->Initialize();
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the UI manager and set verbosities
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // batch mode
  if (!ui) {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro_file);
  } else {
    UImanager->ApplyCommand("/control/execute vis.mac");
    //if ( ui->IsGUI()) {
    //  UImanager->ApplyCommand("/control/execute gui.mac");
   // }
    ui->SessionStart();
    delete ui;
  }

  // stop the timer
  Timer.Stop();
  G4cout << G4endl;
  G4cout << "******************************************";
  G4cout << G4endl;
  G4cout << "Total Real Elapsed Time is: " << Timer.GetRealElapsed();
  G4cout << G4endl;
  G4cout << "Total System Elapsed Time: " << Timer.GetSystemElapsed();
  G4cout << G4endl;
  G4cout << "Total GetUserElapsed Time: " << Timer.GetUserElapsed();
  G4cout << G4endl;
  G4cout << "******************************************";
  G4cout << G4endl;

  delete runManager;
  delete visManager;
  delete workflow_data;

  return 0;
}
