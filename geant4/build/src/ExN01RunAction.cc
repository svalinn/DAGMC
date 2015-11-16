
#include "ExN01RunAction.hh"
#include "ExN01Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01RunAction::ExN01RunAction(UWUW *uwuw_workflow_data)
  : G4UserRunAction()
{
  workflow_data = uwuw_workflow_data;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01RunAction::~ExN01RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file  - DagGeant.root
  G4String fileName = "DagGeant";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01RunAction::EndOfRunAction(const G4Run* run)
{
  // print histogram statistics
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  /*
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
      // convert vol id into string
      std::stringstream ss;
      ss << (it->second).entity_id;
      // get the tally name
      std::string tally_name = ss.str()+"/"+(it->second).tally_type+"/"+(it->second).particle_name;

      int index = 1 + std::distance(workflow_data.tally_library.begin(),it);

      std::cout << index << std::endl;

      // Print out tracklength result per cm2
      G4cout << tally_name+"  : mean = " << analysisManager->GetH1(index)->mean()*cm*cm
                            << " rms = " << analysisManager->GetH1(index)->rms()*cm*cm  << G4endl;

    }
  }
  */
  // need to loop over the histograms and normalise!
  //Double_t norm = hist->GetEntries();
  //hist->Scale(1/norm);
  /*
  for ( it = workflow_data.tally_library.begin() ; it != workflow_data.tally_library.end() ; ++it )
  {
    std::stringstream ss;
    ss << (it->second).entity_id;
    // get the tally name
    std::string tally_name = ss.str()+"_"+(it->second).tally_type+"_"+(it->second).particle_name;
    std::cout << tally_name << std::endl;
    //  G4SDManager::GetSDMpointer()->
  }
  */
  // number of primaries
  //int num_of_event = G4RunManager::GetRunManager()->GetNumberOfEvent();
  int num_of_event = run->GetNumberOfEvent();
  // iterate over tallies
  std::map<std::string,pyne::Tally>::iterator it;
  // loop over histograms and get data
  for ( it = workflow_data->tally_library.begin() ; it != workflow_data->tally_library.end() ; ++it ) {
    int index = 1 + std::distance(workflow_data->tally_library.begin(),it);

    // loop over this histograms
    //    G4cout << index << G4endl;
    G4cout << index << " " << analysisManager->GetH1(index)->mean() << G4endl;
    // scale the result by 1/volume and by (1/cm*cm)
    //analysisManager->ScaleH1(index,1./ double(num_of_event));
  }
  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();

}
