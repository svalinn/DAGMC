
#include "ExN01RunAction.hh"
#include "ExN01Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern UWUW workflow_data;

ExN01RunAction::ExN01RunAction(UWUW data)
 : G4UserRunAction()
{ 
  // get uwuw workflow data
  //workflow_data = data;
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in ExNO1Analysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);

  // iterate over tallies
  std::map<std::string,pyne::Tally>::iterator it;

  // create base tuple
  analysisManager->CreateNtuple("DagGeant","TrackL");
  // Creating histograms
  std::cout << "tallies ! yay" << std::endl;
  for ( it = workflow_data.tally_library.begin() ; it != workflow_data.tally_library.end() ; ++it )
    {
      std::cout << (it->first) << std::endl;
      // convert volid into string
      std::stringstream ss;
      ss << (it->second).entity_id;

      std::string tally_name = ss.str()+"/"+(it->second).tally_type+"/"+(it->second).particle_name;
      // create historgram
      analysisManager->CreateH1(tally_name,(it->second).tally_name,100,0.,1.*m);
      // create tuple
      analysisManager->CreateNtupleDColumn(tally_name);
    }
  // finish tuple
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01RunAction::~ExN01RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "DagGeant";
  analysisManager->OpenFile(fileName);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << "\n ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run \n" << G4endl; 
    }
    else {
      G4cout << "for the local thread \n" << G4endl; 
    }
    
    // iterate over tallies
    std::map<std::string,pyne::Tally>::iterator it;

    // loop over histograms and get data
   for ( it = workflow_data.tally_library.begin() ; it != workflow_data.tally_library.end() ; ++it )
    {
      // convert volid into string
      std::stringstream ss;
      ss << (it->second).entity_id;

      std::string tally_name = ss.str()+"/"+(it->second).tally_type+"/"+(it->second).particle_name;

      int index = 1 + std::distance(workflow_data.tally_library.begin(),it);
  
      std::cout << index << std::endl;

      /*
      G4cout << tally_name+"  : mean = " << G4BestUnit(analysisManager->GetH1(index)->mean(), "Length") 
	                    << " rms = " << G4BestUnit(analysisManager->GetH1(index)->rms(),  "Length") << G4endl;
      */

      G4cout << tally_name+"  : mean = " << analysisManager->GetH1(index)->mean()/(mm*mm)
	     << " rms = " << analysisManager->GetH1(index)->rms()/(mm*mm)  << G4endl;

    }
  }    
  // save histograms & ntuple
  //

  analysisManager->Write();
  analysisManager->CloseFile();  

}
