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

#include "uwuw.hpp"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern UWUW workflow_data;

// constructor for event action
// should populate the ids of the sensitive detectors
ExN01EventAction::ExN01EventAction()
  : G4UserEventAction()
{
  // initialise the list of sensitive volume tally id's
  // so that when we get the first hit we can collect the actual
  // values
  for ( int i = 0 ; i < workflow_data.tally_library.size() ; i++)
  {
    tally_ids.push_back(-1);
  }
  //std::map<std::string,pyne::Tally>::iterator it;
  //for( it = workflow_data.tally_library.begin() ; it != workflow_data.tally_library.end() ; it++ )
//    {
      //vol_tl_ids.push_back(-1);
    //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// destructor, do nothing
ExN01EventAction::~ExN01EventAction()
{}

std::vector<int>
ExN01EventAction::GetTallyID()
{
  std::vector<int> sensitive_detector_ids;

  std::map<std::string,pyne::Tally>::iterator it;
  for( it = workflow_data.tally_library.begin() ; it != workflow_data.tally_library.end() ; ++it )
     {
        // name has to match that in detector construction
//        int idx = std::distance(workflow_data.tally_library.begin(),it);

        // convert int to string
        std::stringstream ss;
        ss <<(it->second).entity_id;
        std::string vol_id = ss.str();
        std::string detector_name;

        if ((it->second).entity_type.find("Volume") != std::string::npos &&
         (it->second).tally_type.find("Flux") != std::string::npos )
          detector_name = "vol_"+vol_id+"_flux/"+(it->second).particle_name+"CellFlux";

        sensitive_detector_ids.push_back(G4SDManager::GetSDMpointer()->GetCollectionID(detector_name));
     }
  return sensitive_detector_ids;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// gets the tally ids
void ExN01EventAction::BeginOfEventAction(const G4Event *)
{
  if(tally_ids[0] == -1)
    tally_ids = GetTallyID();
}

// collect up events and do work
void ExN01EventAction::EndOfEventAction(const G4Event *event)
{
  // get the singleton instance of the analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4cout << ">>> Summary of Event " << event->GetEventID() << G4endl;
  // if the ids aren't initialised then return
  if(std::find(tally_ids.begin(), tally_ids.end(),-1) != tally_ids.end())
    return;

  // get the hit collection
  G4HCofThisEvent *HCE = event->GetHCofThisEvent();
  // Hit Collection for the volume
  G4VHitsCollection *HC;
  // loop over each detector
  std::vector<int>::iterator it;
  for ( it = tally_ids.begin() ; it != tally_ids.end() ; ++it )
  {
     HC = 0;
    // get the hit collection for the tally
     HC = HCE->GetHC(*it);
     // if there are hits
     if(HC->GetSize())
      {
        int n_hit = HC->GetSize();
      }
      // fill the histograms with results
      analysisManager->FillH1(*it+1,1.0);
      analysisManager->FillNtupleDColumn(*it,1.0);
  }
  analysisManager->AddNtupleRow();
}
