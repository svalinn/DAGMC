#include "ExN01DetectorMessenger.hh"
#include "ExN01DetectorConstruction.hh"

ExN01DetectorMessenger::ExN01DetectorMessenger(ExN01DetectorConstruction* Det) : G4UImessenger
(), fDetector(Det)
{
  fDirectory = new G4UIdirectory("/det/");
  fDirectory->SetGuidance("UI command to control B field commands");
  // command for the file based field
  fBFieldCmd = new G4UIcmdWithAString("/det/setBFieldFile", this);
  fBFieldCmd->SetGuidance("Set the file to read the B field from");
  fBFieldCmd->SetParameterName("bFieldFile", false);
  fBFieldCmd->AvailableForStates(G4State_Idle);
  fBFieldCmd->SetToBeBroadcasted(false);
  // command for the wire field current
  fWireFieldCurrentCmd = new G4UIcmdWithADoubleAndUnit("/det/setWireCurrent",this);
  fWireFieldCurrentCmd->SetGuidance("Set the paramters for the wire B field");
  fWireFieldCurrentCmd->SetParameterName("WireCurrent",false);
  fWireFieldCurrentCmd->AvailableForStates(G4State_Idle);
  fWireFieldCurrentCmd->SetToBeBroadcasted(false);
  // command for the wire field radius
  fWireFieldRadiusCmd = new G4UIcmdWithADoubleAndUnit("/det/setWireRadius",this);
  fWireFieldRadiusCmd->SetGuidance("Set the radius for the wire B field");
  fWireFieldRadiusCmd->SetParameterName("WireRadius",false);
  fWireFieldRadiusCmd->AvailableForStates(G4State_Idle);
  fWireFieldRadiusCmd->SetToBeBroadcasted(false);
  // command for the wire permeability
  fWireFieldPermCmd = new G4UIcmdWithADoubleAndUnit("/det/setWireMu",this);
  fWireFieldPermCmd->SetGuidance("Set the permeability for the wire");
  fWireFieldPermCmd->SetParameterName("WirePermeability",false);
  fWireFieldPermCmd->AvailableForStates(G4State_Idle);
  fWireFieldPermCmd->SetToBeBroadcasted(false);
  // contruct method - method to tell Geant4 to construct the class
  fFieldConstructCmd = new G4UIcmdWithoutParameter("/det/construct",this);
  fFieldConstructCmd->AvailableForStates(G4State_Idle);
  fFieldConstructCmd->SetToBeBroadcasted(false);

}

ExN01DetectorMessenger::~ExN01DetectorMessenger() {
  delete fBFieldCmd;
  delete fWireFieldCurrentCmd;
  delete fWireFieldRadiusCmd;
  delete fWireFieldPermCmd;
  delete fFieldConstructCmd;
  delete fDirectory;  
}

void ExN01DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  // set the b field filename
  if (command == fBFieldCmd) {
    fDetector->SetBFieldFileName(newValue);
  }
  // we have recieved the wirefield current cmd
  if (command == fWireFieldCurrentCmd) {
    G4double value = G4UIcmdWithADoubleAndUnit::GetNewDoubleRawValue(newValue);
    if ( value == 0.0 ) {
      G4cout << "Current set to 0, must be +ve or -ve" << G4endl; 
    }
    fDetector->SetWireFieldCurrent(value);
  }
  // we have recieved the radius command
  if (command == fWireFieldRadiusCmd) {
    G4double value = G4UIcmdWithADoubleAndUnit::GetNewDoubleRawValue(newValue);
    if ( value <= 0.0 ) {
      G4cout << "Wire Radius 0 or less, must be +ve" << G4endl; 
    }
    fDetector->SetWireFieldRadius(G4UIcmdWithADoubleAndUnit::GetNewDoubleRawValue(newValue));
  }
  // we have recieved the permeability command
  if (command == fWireFieldPermCmd) {
    G4double value = G4UIcmdWithADoubleAndUnit::GetNewDoubleRawValue(newValue);
    if ( value <= 0.0 ) {
      G4cout << "Wire permeability 0 or less, must be +ve" << G4endl; 
    }    
    fDetector->SetWireFieldPermeability(G4UIcmdWithADoubleAndUnit::GetNewDoubleRawValue(newValue));
  }
  // we have recieved the contruct command
  if (command == fFieldConstructCmd) {
    fDetector->ConstructSDandField();
  }
}

G4String ExN01DetectorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String ans;
  if (command == fBFieldCmd) {
    ans = fDetector->GetBFieldFileName();
  }
  if (command == fWireFieldCurrentCmd) {
    ans = fDetector->GetWireFieldCurrent();
  }
  if (command == fWireFieldRadiusCmd) {
    ans = fDetector->GetWireFieldRadius();
  }
  if (command == fWireFieldPermCmd) {
    ans = fDetector->GetWireFieldPerm();
  }
  return ans;
}
