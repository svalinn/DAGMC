#ifndef ExN01DetectorMessenger_hh
#define ExN01DetectorMessenger_hh 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"

#include "globals.hh"

class ExN01DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class ExN01DetectorMessenger : public G4UImessenger
{
  public:
    ExN01DetectorMessenger(ExN01DetectorConstruction*);
    ~ExN01DetectorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);
    virtual G4String GetCurrentValue(G4UIcommand* command);

  private:
    ExN01DetectorConstruction* fDetector;

    G4UIdirectory* fDirectory;
    G4UIcmdWithoutParameters* fFieldConstructCmd;
    G4UIcmdWithAString* fBFieldCmd;
    G4UIcmdWithADoubleAndUnit* fWireFieldCurrentCmd;
    G4UIcmdWithADoubleAndUnit* fWireFieldRadiusCmd;
    G4UIcmdWithADoubleAndUnit* fWireFieldPermCmd;
  
};

#endif
