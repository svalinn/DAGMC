#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "ExN01Run.hh"

ExN01Run::ExN01Run() : G4Run() {  
    std::vector<G4String> detector_types = {"Dose","Flux","NoSecondary","NoStep"};
    std::vector<G4String> particle_names = {"Neutron","Photon","Electron"};
    
    //G4String det_name = "Dose"; // maybe add more in the future?
    // get sensitive detector manager
    G4SDManager *SDM = G4SDManager::GetSDMpointer();
    G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
    std::vector<G4LogicalVolume*>::iterator lv_it;
    std::vector<G4String>::iterator det_it;
    std::vector<G4String>::iterator par_it;
    // loop over all the logical volumes
    for ( lv_it = lvs->begin(); lv_it != lvs->end(); lv_it++ ) {
      for ( det_it = detector_types.begin() ; det_it != detector_types.end() ; det_it++ ) {
	for ( par_it = particle_names.begin() ; par_it != particle_names.end() ;  par_it++ ) {
	  G4String lv_name = (*lv_it)->GetName();
	  // make the name
	  G4String det_name = *det_it; // name of score type
	  G4String particle_name = *par_it; // name of particle
	  G4String detector_name = lv_name + "/" + det_name + "/" + particle_name;   
	  // get unique collection id 
	  G4int score_id = SDM->GetCollectionID(detector_name);    
	  if ( score_id != -1 ) {
	    detectors[detector_name] = score_id;
	  }
	}
      }	
    }
    // all done
}

ExN01Run::~ExN01Run(){}

void ExN01Run::RecordEvent(const G4Event *evt) {
  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(!HCE) return;
  numberOfEvent++;
  // accumlate scores for this event
  std::map<G4String,G4int>::iterator it;
  for ( it = detectors.begin() ; it != detectors.end() ; ++it ) {
    // unqiue score id
    G4int score_id = it->second;
    // the actual tally score
    G4THitsMap<G4double>* score = (G4THitsMap<G4double>*)HCE->GetHC(score_id);
    // score sum 
    total_score[score_id] += *score;
    //G4THitsMap<G4double> val += *score;
  }  
}

void ExN01Run::Merge(const G4Run *aRun) {
  G4Run::Merge(aRun);
}

/*
// Get the score for the given id
G4THitsMap<G4double> ExN01Run::GetScore(G4int id) {
  G4THitsMap<G4double> tot = total_score[id];
  return tot;
}

// get the score for a given arbitrary volume and score name
G4THitsMap<G4double> ExN01Run::GetScore(G4LogicalVolume *vol, G4String score_name) {
  G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
  G4SDManager *SDM = G4SDManager::GetSDMpointer();

  std::vector<G4LogicalVolume*>::iterator lv_it;
  // loop over all the logical volumes
  for ( lv_it = lvs->begin(); lv_it != lvs->end(); ++lv_it ) {
    if( vol == *lv_it) { 
      G4String lv_name = (*lv_it)->GetName();
      // make the name
      // G4String det_name = "Dose";
      G4String detector_name = lv_name + "/" + score_name;
      // get unique collection id 
      G4int score_id = SDM->GetCollectionID(detector_name);
      if(score_id != -1 ) 
        return GetScore(score_id);
      else {
        G4THitsMap<G4double> dummy;
        return dummy;
      }
    }
  }
}

// Get total by score id
G4double ExN01Run::GetTotal(G4int id) {
  G4THitsMap<G4double> scores = GetScore(id);
  G4double total = 0.;
  if(scores.GetSize()==0) return total;
  std::map<G4int,G4double*>::iterator itr = scores.GetMap()->begin();
  for(; itr != scores.GetMap()->end(); itr++) { 
    total += *(itr->second); 
  }
  total /= double(numberOfEvent);
  return total;
}
*/

// Get total by logical volume and score name
G4double ExN01Run::GetTotal(G4LogicalVolume *vol, G4String score_name) {

  G4String lv_name = vol->GetName();
  // make the name
  // G4String det_name = "Dose";
  G4String detector_name = lv_name + "/" + score_name;
  // get unique collection id 
  G4SDManager *SDM = G4SDManager::GetSDMpointer();
  G4int score_id = SDM->GetCollectionID(detector_name);

  if ( score_id < 0 ) return 0.0;
  G4double total = 0.;
  //G4cout << " $" << score_id << " " << total_score[score_id].GetSize() << G4endl;

  //  G4THitsMap<G4double> scores = total_score[score_id];
  if(total_score[score_id].GetSize() == 0 ) return total;
  //  return total;
  std::map<G4int,G4double*>::iterator it;
  for ( it = total_score[score_id].GetMap()->begin() ; it != total_score[score_id].GetMap()->end() ; it++ ) {
    total += *(it->second);
  }
  /* note very well that in the above loop, it would be rational to define a 
   * temporary variable like
   *  G4THitsMap<G4double> scores = total_score[score_id];
   * if this is done and the iterator loops over scores.GetMap() this leads to a
   * segfault at program termination. Maybe this variable is stored on the stack 
   * rather than the heap. ? Don't really know, will report to Geant4 team.
   */
  total /= double(numberOfEvent);
  return total;
}
