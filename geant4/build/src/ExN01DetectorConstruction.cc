//
//
// $Id$
//

#include "ExN01DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"

#include "DagSolid.hh"
#include "MBInterface.hpp"
#include "DagMC.hpp"

#include "DagSolidMaterial.hh"
#include "DagSolidTally.hh"

using namespace moab;

 // dag_volumes
DagMC* dagmc = DagMC::instance(); // create dag instance
 

ExN01DetectorConstruction::ExN01DetectorConstruction(std::string uwuw_file)
  :  world_volume_log(0)

{
  uwuw_filename = uwuw_file;
  ;
}

ExN01DetectorConstruction::~ExN01DetectorConstruction()
{
}

G4VPhysicalVolume* ExN01DetectorConstruction::Construct()
{

  // const char* h5mfilename = "/data/opt/DagGeant4/atic_uwuw_zip.h5m";
 
  //  std::string dag_file(uwuw_filename);


  // check for Materials first
  try
    {
      pyne::Material mat;
      mat.from_hdf5(uwuw_filename,"/materials");
    }
  catch (const std::exception &except) // catch the exception from from_hdf5                                                      
    {
      std::cout << "No Materials found in the file, " << uwuw_filename << std::endl;
      std::cout << "Cannot continue without materials " << std::endl;

      exit(1);
    }


  // load the material from the UW^2 library
  std::map<std::string,G4Material*> material_lib;
  //  material_lib = load_uwuw_materials(dag_file);

  material_lib = load_uwuw_materials(uwuw_filename);

  //  std::cout << material_lib["mat:\Lead"] << std::endl;
    
  G4VisAttributes * invis = new G4VisAttributes(G4VisAttributes::Invisible);

  //------------------------------------------------------ volumes
  // -- World Volume in which we place other volumes
  
  
  G4double world_width  = 1020.0*cm;
  G4double world_height = 1020.0*cm;
  G4double world_depth  = 1020.0*cm;
  
  G4Box* world_volume = new G4Box("world_volume_box",world_width,world_height,world_depth);
  world_volume_log = new G4LogicalVolume(world_volume,material_lib["mat:Vacuum"],"world_vol_log",0,0,0);
  world_volume_log->SetVisAttributes(invis);
  G4PVPlacement* world_volume_phys = new G4PVPlacement(0,G4ThreeVector(),world_volume_log,
						      "world_vol",0,false,0); 
   
  G4cout << "Load sample file = "     << dagmc->load_file(uwuw_filename.c_str(),0) << G4endl;
  G4cout << "Initialize OBB tree = "  << dagmc->init_OBBTree() << G4endl;

  G4int Num_of_survertices = dagmc->num_entities(2);
  G4int num_of_objects = dagmc->num_entities(3);

  G4cout << "There are " << num_of_objects << " dag volumes" << G4endl;

  // parse the dagmc props
  std::vector<std::string> dagsolid_keywords;
  dagsolid_keywords.push_back("mat");
  dagsolid_keywords.push_back("rho");
  std::map<std::string,std::string> synonymns;
  char *delimeters = ":/";

  MBErrorCode rval = dagmc->parse_properties(dagsolid_keywords,synonymns,delimeters);
  
  //Store a list of DagSolids, Logical Vols, and Physical Vols
  std::vector<DagSolid*> dag_volumes;
  std::vector<G4LogicalVolume*> dag_logical_volumes;
  std::vector<G4PVPlacement*> dag_physical_volumes;
  
  for ( int dag_idx = 1 ; dag_idx < num_of_objects ; dag_idx++ )
    {
      // get the MBEntity handle for the volume
      int dag_id = dagmc->id_by_index(3,dag_idx);
      MBEntityHandle volume = dagmc->entity_by_id(3,dag_id);
      std::string dag_material_name;
      rval = dagmc->prop_value(volume,"mat",dag_material_name);
      std::cout << "id = " << dag_id << " mat = " << dag_material_name << std::endl;
      std::stringstream int_to_string;
      int_to_string << dag_idx;
      std::string idx_str = int_to_string.str();

      // create new volume
      DagSolid* dag_vol = new DagSolid("vol_"+idx_str,dagmc,dag_idx);
      dag_volumes.push_back(dag_vol);

      // define logical volume
      if( dag_material_name.find("Blackhole") != std::string::npos || 
	  dag_material_name.find("Graveyard") != std::string::npos )
	dag_material_name = "Vacuum";

      G4LogicalVolume* dag_vol_log = new G4LogicalVolume(dag_vol,material_lib["mat:"+dag_material_name],
							"vol_"+idx_str+"_log",0,0,0);
      
      std::cout << "vol_" << idx_str << "  has prop " << material_lib["mat:"+dag_material_name] << std::endl;
      dag_logical_volumes.push_back(dag_vol_log);

      // define physical volumes
      G4PVPlacement* dag_vol_phys = new G4PVPlacement(0,G4ThreeVector(0*cm,0*cm,0*cm),dag_vol_log,
						     "volume_"+idx_str+"_phys",
						      world_volume_log,false,0);
      dag_physical_volumes.push_back(dag_vol_phys);
    }
  

  return world_volume_phys;
}


void ExN01DetectorConstruction::ConstructSDandField()
{
  
  //  std::string filename = "atic_uwuw_zip.h5m";
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // check for tallies first
  try
    {
      pyne::Tally tally;
      tally.from_hdf5(uwuw_filename,"/tally");
    }
  catch (const std::exception &except) // catch the exception from from_hdf5                                                      
    {
      std::cout << "No Tallies found in the file, " << uwuw_filename << std::endl;
      return;
    }

  // load 
  std::map<std::string,pyne::Tally>::iterator t_it;

  //  load tallies from the h5m file
  tally_library = load_uwuw_tallies(uwuw_filename);
  
  G4VPrimitiveScorer* primitive; //set g4scorer

  // TrackLength Scorer   
  std::vector<G4MultiFunctionalDetector*> detectors;

  // loop over the tallies
  for ( t_it = tally_library.begin() ; t_it != tally_library.end() ; ++t_it ) 
    {
      pyne::Tally scorer = t_it -> second; // current tally
      if( scorer.tally_type.find("Flux") != std::string::npos && 
	  scorer.entity_type.find("Volume") != std::string::npos )
	{
	  int vol_id = scorer.entity_id;
	  MBEntityHandle vol = dagmc->entity_by_id(3,vol_id); // convert id to eh
	  int vol_idx = dagmc->index_by_handle(vol); // get the volume index 
	  // now get logical volume ( vol_idx - 1 )
	  
	  // convert vol_idx to string
	  std::stringstream int_to_string;
	  int_to_string << vol_idx;
	  std::string idx_str = int_to_string.str();
	  
	  // create G4Multi
	  G4MultiFunctionalDetector* vol_det = new G4MultiFunctionalDetector("vol_"+idx_str+"_flux");
	  primitive = new G4PSTrackLength("TrackLength");
	  vol_det->RegisterPrimitive(primitive);
	  SetSensitiveDetector("vol_"+idx_str+"_log",vol_det);
      
	  // store detectors
	  detectors.push_back(vol_det);
	}
    }
  
}
