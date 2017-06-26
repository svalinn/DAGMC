#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "ExN01Run.hh"

ExN01Run::ExN01Run() : nEvent(0) {
    // get sensitive detector manager
    G4SDManager *SDM = G4SDManager::GetSDMpointer();
    
    G4String det_name = "Dose"; // maybe add more in the future?

    G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
    std::vector<G4LogicalVolume*>::iterator lv_it;
    // loop over all the logical volumes
    for ( lv_it = lvs->begin(); lv_it != lvs->end(); lv_it++ ) {
      G4String lv_name = (*lv_it)->GetName();
      // make the name
      G4String detector_name = lv_name + "/" + det_name;
      // get unique collection id 
      G4int score_id = SDM->GetCollectionID(detector_name);
      detectors[detector_name] = score_id;
    }
    // all done
}

ExN01Run::~ExN01Run(){}

void ExN01Run::RecordEvent(const G4Event *evt) {
  nEvent++;
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  // accumlate scores for this event
  std::map<G4String,G4int>::iterator it;
  for ( it = detectors.begin() ; it != detectors.end() ; it++ ) {
    // unqiue score id
    G4int score_id = it->second;
    // the actual tally score
    G4THitsMap<G4double>* score = (G4THitsMap<G4double>*)HCE->GetHC(score_id);
    // score sum 
    total_score[it->second] += *score;
  }  
}

void ExN01Run::Merge(const G4Run *aRun) {
  G4Run::Merge(aRun);
}

G4THitsMap<G4double> ExN01Run::GetTotal(G4int id) {
  G4THitsMap<G4double> tot = total_score[id];
  return tot;
}

G4THitsMap<G4double> ExN01Run::GetTotal(G4LogicalVolume *vol, G4String score_name) {
  G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
  G4SDManager *SDM = G4SDManager::GetSDMpointer();

  std::vector<G4LogicalVolume*>::iterator lv_it;
  // loop over all the logical volumes
  for ( lv_it = lvs->begin(); lv_it != lvs->end(); lv_it++ ) {
    if( vol == *lv_it) { 
      G4String lv_name = (*lv_it)->GetName();
      // make the name
      G4String det_name = "Dose";
      G4String detector_name = lv_name + "/" + det_name;
      // get unique collection id 
      G4int score_id = SDM->GetCollectionID(detector_name);
      // 
      return GetTotal(score_id);
    }
  }
//  return G4THitsMap<G4double>*;
}
