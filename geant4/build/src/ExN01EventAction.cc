//
/// ile ExN01EventAction.cc
/// rief Implementation of the ExN01EventAction class

#include "ExN01EventAction.hh"
#include "ExN01Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "../pyne/pyne.h"
#include "ExN01DetectorConstruction.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// constructor for event action
// should populate the ids of the sensitive detectors
ExN01EventAction::ExN01EventAction()
  : G4UserEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// destructor, do nothing
ExN01EventAction::~ExN01EventAction()
{}

// gets the tally ids
void ExN01EventAction::BeginOfEventAction(const G4Event *)
{
}

// collect up events and do work
void ExN01EventAction::EndOfEventAction(const G4Event *event)
{
  // get the stored trajectories
  /*
  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  // periodic printing
  G4int eventID = event->GetEventID();
  if ( eventID < 100 || eventID % 100 == 0 ) {
    G4cout << ">>> Event: " << eventID << G4endl;
    if ( trajectoryContainer ) {
      G4cout << "    " << n_trajectories
             << " trajectories stored in this event" << G4endl;
    }
    G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
    G4cout << "     "
           << hc->GetSize() << " hits stored in this event" << G4endl;
  }
  */
  /* will need to setup histograms
  // get the singleton instance of the analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

       // fill the histograms with results
       analysisManager->FillH1(*it+1,1.0);
       analysisManager->FillNtupleDColumn(*it,1.0);
   }
   analysisManager->AddNtupleRow();
   */
}
