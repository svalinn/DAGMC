
#include "ExN01RunAction.hh"

#include "ExN01Analysis.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

ExN01RunAction::ExN01RunAction(UWUW* uwuw_workflow_data) : G4UserRunAction() {
  workflow_data = uwuw_workflow_data;
  // create an analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
}

ExN01RunAction::~ExN01RunAction() {}


void ExN01RunAction::BeginOfRunAction(const G4Run* /*run*/) {
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();


  // Open an output file
  // G4String fileName = "DagGeant";
  analysisManager->OpenFile();
  return;
}
void ExN01RunAction::EndOfRunAction(const G4Run* run) {
  // print histogram statistics
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Save histograms
  analysisManager->Write();
  analysisManager->CloseFile();
  return;

}
