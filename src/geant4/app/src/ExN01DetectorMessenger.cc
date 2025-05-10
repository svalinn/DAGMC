#include "ExN01DetectorMessenger.hh"
#include "ExN01DetectorConstruction.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"

ExN01DetectorMessenger::ExN01DetectorMessenger(ExN01DetectorConstruction* Det) : G4UImessenger
(), fDetector(Det)
{
  fDirectory = new G4UIdirectory("/det/");
  fDirectory->SetGuidance("UI command to control B field commands");
  fBFieldCmd = new G4UIcmdWithAString("/det/setBFieldFile", this);
  fBFieldCmd->SetGuidance("Set the file to read the B field from");
  fBFieldCmd->SetParameterName("bFieldFile", false);
  fBFieldCmd->AvailableForStates(G4State_Idle);
  fBFieldCmd->SetToBeBroadcasted(false);
}

ExN01DetectorMessenger::~ExN01DetectorMessenger() {
  delete fBFieldCmd;
  delete fDirectory;  
}

void ExN01DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == fBFieldCmd) {
    fDetector->SetBFieldFileName(newValue);
  }
}

G4String ExN01DetectorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String ans;
  if (command == fBFieldCmd) {
    ans = fDetector->GetBFieldFileName();
  }
  return ans;
}