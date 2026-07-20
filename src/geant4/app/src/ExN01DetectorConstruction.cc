//
//
// $Id$
//

#include "ExN01DetectorConstruction.hh"

#include "DagSolid.hh"
#include "DagSolidColors.hh"
#include "DagSolidMaterial.hh"
#include "DagSolidTally.hh"
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

G4ThreadLocal MagneticField* ExN01DetectorConstruction::fMagneticField = nullptr;
G4ThreadLocal G4FieldManager* ExN01DetectorConstruction::fFieldMgr = nullptr;

// constructor
ExN01DetectorConstruction::ExN01DetectorConstruction(UWUW* uwuw_workflow_data, 
                                                     G4String magnetic_field_filename)
  : fWorldVolumeLog(0),fDetectorMessenger(0),fWireCurrent(0),
    fWireRadius(0),fWirePermeability(0) {
  workflow_data = uwuw_workflow_data;
  fBfieldFilename = magnetic_field_filename;

  dagmc = new moab::DagMC();
  // load the material from the UW^2 library
  material_lib = load_uwuw_materials(workflow_data);
  fDetectorMessenger = new ExN01DetectorMessenger(this);
}

// destructor
ExN01DetectorConstruction::~ExN01DetectorConstruction() {}

// the main method - takes the problem and loads
G4VPhysicalVolume* ExN01DetectorConstruction::Construct() {

#ifdef GEANT4_GT_11
  G4VisAttributes* invis = new G4VisAttributes(G4VisAttributes::GetInvisible());
#else
  G4VisAttributes* invis = new G4VisAttributes(G4VisAttributes::Invisible);
#endif

  //------------------------------------------------------ volumes
  // -- World Volume in which we place other volumes
  // we can probably do something like get the obb for the implicit complement
  // and use that as the size for the world vol

  // load the dagmc file - only to get the counts of volumes
  moab::ErrorCode rval = dagmc->load_file(workflow_data->full_filepath.c_str());
  if (rval != moab::MB_SUCCESS) {
    G4cout << "ERROR: Failed to load the DAGMC file " << workflow_data->full_filepath << G4endl;
    exit(1);
  }

  // build the trees
  rval = dagmc->init_OBBTree();
  if (rval != moab::MB_SUCCESS) {
    G4cout << "ERROR: Failed to build the OBB Tree" << G4endl;
    exit(1);
  }

  // attach a metadata instance
  DMD = new dagmcMetaData(dagmc);
  DMD->load_property_data();

  // set the world
  G4double world_width = GetMaxOrdinate() * cm;
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2. * world_width);

  G4Box* world_volume =
      new G4Box("world_volume_box", world_width, world_width, world_width);
  fWorldVolumeLog = new G4LogicalVolume(
      world_volume, material_lib["mat:Vacuum"], "world_vol_log", 0, 0, 0);
  fWorldVolumeLog->SetVisAttributes(invis);
  G4PVPlacement* world_volume_phys = new G4PVPlacement(
      0, G4ThreeVector(), fWorldVolumeLog, "world_vol", 0, false, 0);

  // get count of entities
  G4int num_of_objects = dagmc->num_entities(3);

  // set the colours to be used
  UniformColorGenerator* colorgen =
      new UniformColorGenerator(material_lib.size());
  // UniformColorGenerator *colorgen = new
  // UniformColorGenerator(num_of_objects);
  colorgen->Generate();
  std::vector<RGB> colours = colorgen->GetColors();

  G4cout << "There are " << num_of_objects << " dag volumes" << G4endl;

  // Store a list of DagSolids, Logical Vols, and Physical Vols
  std::vector<DagSolid*> dag_volumes;
  std::vector<G4PVPlacement*> dag_physical_volumes;

  // load the properties from the metadata instance
  for (int dag_idx = 1; dag_idx < num_of_objects; dag_idx++) {
    G4String idx_str = std::to_string(dag_idx);
    // get the MBEntity handle for the volume
    int dag_id = dagmc->id_by_index(3, dag_idx);
    moab::EntityHandle volume = dagmc->entity_by_id(3, dag_id);
    // get the material_name
    std::string mat_name = DMD->volume_material_property_data_eh[volume];

    // create new volume
    DagSolid* dag_vol = new DagSolid("vol_" + idx_str, dagmc, dag_idx);
    dag_volumes.push_back(dag_vol);
    // make new logical volume
    std::string material_name = mat_name;
    if (mat_name == "mat:Graveyard" || mat_name == "mat:Vacuum") {
      material_name = "mat:Vacuum";
    }

    G4LogicalVolume* dag_vol_log =
        new G4LogicalVolume(dag_vol, material_lib[material_name],
                            "vol_" + idx_str + "_log", 0, 0, 0);
    // set the vis attributes
    if (mat_name == "mat:Graveyard" || mat_name == "mat:Vacuum") {
      dag_vol_log->SetVisAttributes(invis);
    } else {
      // use distance in place of an index
      G4int mat_idx =
          std::distance(material_lib.begin(), material_lib.find(mat_name));
      // set the colour
      dag_vol_log->SetVisAttributes(
          G4Color(colours[mat_idx].r, colours[mat_idx].g, colours[mat_idx].b));
    }
    dag_logical_volumes.push_back(dag_vol_log);
    // make a new physical placement
    G4PVPlacement* dag_vol_phys = new G4PVPlacement(
        0, G4ThreeVector(0 * cm, 0 * cm, 0 * cm), dag_vol_log,
        "volume_" + idx_str + "_phys", fWorldVolumeLog, false, 0);
    dag_physical_volumes.push_back(dag_vol_phys);
  }

  return world_volume_phys;
}

// find the maximum coordinate in any direction using the graveyard
// used to set a sensible mother volume size
G4double ExN01DetectorConstruction::GetMaxOrdinate() {
  G4double max_ordinate = 0.0;  
  // find the max ordinate by looping over each volume
  for (int dag_idx = 1; dag_idx < dagmc->num_entities(3); dag_idx++) {
    int dag_id = dagmc->id_by_index(3, dag_idx);
    moab::EntityHandle volume = dagmc->entity_by_id(3, dag_id);
    // get the material_name
    G4double min[3], max[3];
    moab::ErrorCode rval = dagmc->getobb(volume, min, max);
    // find the max ordinate
    G4double x = std::max(std::abs(min[0]), std::abs(max[0]));
    G4double y = std::max(std::abs(min[1]), std::abs(max[1]));
    G4double z = std::max(std::abs(min[2]), std::abs(max[2]));
  
    G4double maxord = std::max(x, std::max(y, z));

    // if max is bigger than any other found
    if (maxord > max_ordinate) {
      max_ordinate = maxord;
    }
  }

  // return a value 20% bigger 
  return max_ordinate*1.2;
}

void ExN01DetectorConstruction::SetBFieldFileName(G4String filename) {
  // set the filename
  fBfieldFilename = filename;
}

G4String ExN01DetectorConstruction::GetBFieldFileName() {
  // set the filename
  return fBfieldFilename;
}

void ExN01DetectorConstruction::SetWireFieldCurrent(const G4double current) {
  // set the current
  fWireCurrent = current;
}

G4double ExN01DetectorConstruction::GetWireFieldCurrent() {
  // set the filename
  return fWireCurrent;
}

void ExN01DetectorConstruction::SetWireFieldRadius(const G4double radius) {
  // set the current
  fWireRadius = radius;
}

G4double ExN01DetectorConstruction::GetWireFieldRadius() {
  // set the filename
  return fWireRadius;
}

void ExN01DetectorConstruction::SetWireFieldPermeability  (const G4double permeability) {
  // set the current
  fWirePermeability = permeability;
}

G4double ExN01DetectorConstruction::GetWireFieldPermeability() {
  // set the filename
  return fWirePermeability;
}

void ExN01DetectorConstruction::ConstructSDandField() {
  // instanciate the magnetic field
  if (!fBfieldFilename.empty()) {
    G4cout << "Loading magnetic field from file: " << fBfieldFilename
           << "..." << G4endl;
    if (fMagneticField) {
      delete fMagneticField;
      fMagneticField = nullptr;
    }
    fMagneticField = new MagneticField();
    fMagneticField->LoadFile(fBfieldFilename);
    if(fFieldMgr) {
      delete fFieldMgr;
      fFieldMgr = nullptr;
    }
    fFieldMgr = new G4FieldManager(fMagneticField);
    fFieldMgr->SetDetectorField(fMagneticField);
    fFieldMgr->CreateChordFinder(fMagneticField);
    // set the field to be global
    fWorldVolumeLog->SetFieldManager(fFieldMgr, true);
    G4cout << "Magnetic field loaded." << G4endl;
  } else {
    G4cout << "No magnetic field file specified." << G4endl;
    // clear out any existing field
    if (fMagneticField) {
      delete fMagneticField;
      fMagneticField = nullptr;
    }

    // check for issues
    if ( fWireCurrent == 0. ) {
      G4Exception("Ex01DetectorConstruction","1",FatalErrorInArgument,"Wire Current is 0"); 
    }
    // make sure that radius is gt 0
    if ( fWireRadius <= 0. ) {
      G4Exception("Ex01DetectorConstruction","1",FatalErrorInArgument,"Wire Radius <= 0");       
    }
    // make sure permeability
    if ( fWirePermeability <= 0. ) { 
      G4Exception("Ex01DetectorConstruction","1",FatalErrorInArgument,"Wire Permeability <= 0");      
    }

    // create a new wire magnetic field
    fMagneticField = new WireMagneticField(fWireCurrent,
					   fWireRadius,
					   fWirePermeability);
    // attach to the field manager
    fFieldMgr = new G4FieldManager(fMagneticField);
    fFieldMgr->SetDetectorField(fMagneticField);
    fFieldMgr->CreateChordFinder(fMagneticField);
    // set the field to be global
    fWorldVolumeLog->SetFieldManager(fFieldMgr, true);
    G4cout << "Using a simple wire based magnetic field." << G4endl;
  }
} 
