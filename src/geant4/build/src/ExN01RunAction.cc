#include "ExN01Run.hh"
#include "ExN01RunAction.hh"
#include "ExN01Analysis.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01RunAction::ExN01RunAction(UWUW* uwuw_workflow_data) : G4UserRunAction() {
  workflow_data = uwuw_workflow_data;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01RunAction::~ExN01RunAction() {
}

G4Run* ExN01RunAction::GenerateRun() { 
  return new ExN01Run(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01RunAction::BeginOfRunAction(const G4Run* /*run*/) {
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file  - DagGeant.root
  G4String fileName = "DagGeant";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01RunAction::EndOfRunAction(const G4Run* run) {
  // print histogram statistics
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // number of primaries
  // cast to ExN01Run
  ExN01Run* theRun = (ExN01Run*) run;
  
  G4cout << "############################################################################" << G4endl;
  G4cout << "#                                                                          #" << G4endl;
  G4cout << "#                           End of Run Summary                             #" << G4endl;
  G4cout << "#                                                                          #" << G4endl;
  G4cout << "############################################################################" << G4endl;
  G4cout << "#                                                                          #" << G4endl;
  G4cout << "#              Volume Name           Score                                 #" << G4endl;

  // loop over the logical volumes and print summary
  G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
  std::vector<G4LogicalVolume*>::iterator lv_it;
  for( lv_it = lvs->begin(); lv_it != lvs->end(); ++lv_it ) {
    G4String lv_name = (*lv_it)->GetName();
    G4cout << "#        " << lv_name;
    G4String score_type = "Dose";
    //G4double score = 0.0;
    G4double score = theRun->GetTotal( *lv_it, score_type);
    G4cout << "      " << score/MeV*g << G4endl;
  }

  G4cout << "#                                                                          #" << G4endl;
  G4cout << "############################################################################" << G4endl;
  
  // iterate over tallies
  std::map<std::string, pyne::Tally>::iterator it;
  // loop over histograms and get data
  for (it = workflow_data->tally_library.begin() ; it != workflow_data->tally_library.end() ; ++it) {
    int index = 1 + std::distance(workflow_data->tally_library.begin(), it);

    // loop over this histograms
    G4cout << index << " " << analysisManager->GetH1(index)->mean() << G4endl;
  }
  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();

}
