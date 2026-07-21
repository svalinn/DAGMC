
#ifndef ExN01DetectorConstruction_HH
#define ExN01DetectorConstruction_HH 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include <map>
#include <string>
#include <vector>

#include "DagMC.hpp"
#include "DagSolid.hh"
#include "DagSolidColors.hh"
#include "DagSolidMagneticField.hh"
#include "DagSolidMaterial.hh"
#include "DagSolidTally.hh"

#include "ExN01DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "dagmcmetadata.hpp"
#include "moab/Interface.hpp"
#include "pyne.h"

#ifndef uwuw_hpp
#define uwuw_hpp 1
#include "uwuw.hpp"
#endif

class ExN01DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  ExN01DetectorConstruction(UWUW* uwuw_workflow_data, G4String magnetic_field_filename);
  ~ExN01DetectorConstruction();

 public:
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  void SetBFieldFileName(G4String filename);
  void SetWireFieldCurrent(G4double current);
  void SetWireFieldRadius(const G4double radius);
  void SetWireFieldPermeability(const G4double permeability);
  G4String GetBFieldFileName();
  G4double GetWireFieldCurrent();
  G4double GetWireFieldRadius();
  G4double GetWireFieldPermeability();

  G4double GetMaxOrdinate();  

 private:
  static G4ThreadLocal G4MagneticField* fMagneticField;
  static G4ThreadLocal G4FieldManager* fFieldMgr;

  // Logical volumes
  //
  //G4MagneticField* magField;
  G4LogicalVolume* fWorldVolumeLog;
  //std::map<int, G4LogicalVolume*> dag_logical_volumes;
  std::vector<G4LogicalVolume*> dag_logical_volumes;
  UWUW* workflow_data;
  moab::DagMC* dagmc;
  dagmcMetaData* DMD;
  std::map<std::string, G4Material*> material_lib;
  G4String fBfieldFilename;
  G4double fWireCurrent;
  G4double fWireRadius;
  G4double fWirePermeability;
  ExN01DetectorMessenger *fDetectorMessenger;
};

#endif
