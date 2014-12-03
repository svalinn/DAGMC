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

#include "pyne/pyne.h"
#include "ExN01DetectorConstruction.hh"

#include "Randomize.hh"
#include <iomanip>

#include "uwuw.hpp"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern UWUW workflow_data;

ExN01EventAction::ExN01EventAction()
  : G4UserEventAction()
{
  // build the tl_id vector
  std::map<std::string,pyne::Tally>::iterator it;
  for( it = workflow_data.tally_library.begin() ; it != workflow_data.tally_library.end() ; it++ ) 
    {
      vol_tl_ids.push_back(-1);
      tracklengths.push_back(0.0);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01EventAction::~ExN01EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4THitsMap<G4double>* 
ExN01EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  G4THitsMap<G4double>* hitsCollection 
    = static_cast<G4THitsMap<G4double>*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("ExN01EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ExN01EventAction::GetSum(G4THitsMap<G4double>* hitsMap) const
{
  G4double sumValue = 0;
  std::map<G4int, G4double*>::iterator it;
  for ( it = hitsMap->GetMap()->begin(); it != hitsMap->GetMap()->end(); it++) {
    sumValue += *(it->second);
  }
  return sumValue;  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01EventAction::EndOfEventAction(const G4Event* event)
{  
  // get analysis manager                                                                                                          
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


  std::map<std::string,pyne::Tally>::iterator it;
  if ( vol_tl_ids[0] == -1 )
    {
      for( it = workflow_data.tally_library.begin() ; it != workflow_data.tally_library.end() ; ++it )
	{
	  // name has to match that in detector construction 
          int idx = std::distance(workflow_data.tally_library.begin(),it);

	  // convert int to string
	  std::stringstream ss;
	  ss <<(it->second).entity_id;

	  std::string vol_id = ss.str();

	  std::string detector_name;
	  if ((it->second).entity_type.find("Volume") != std::string::npos && 
	      (it->second).tally_type.find("Flux") != std::string::npos )
	    detector_name = "vol_"+vol_id+"_flux/"+(it->second).particle_name+"CellFlux";

	  vol_tl_ids[idx] = G4SDManager::GetSDMpointer()->GetCollectionID(detector_name);
	}
    }
  
  for( it = workflow_data.tally_library.begin() ; it != workflow_data.tally_library.end() ; ++it )
    {	  
      int idx = std::distance(workflow_data.tally_library.begin(),it);
      tracklengths[idx] = GetSum(GetHitsCollection(vol_tl_ids[idx], event));
      analysisManager->FillH1(idx+1, tracklengths[idx]);
      analysisManager->FillNtupleDColumn(idx,tracklengths[idx]);
    }
  analysisManager->AddNtupleRow();
}
  

