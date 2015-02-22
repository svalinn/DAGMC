//
//
// $Id$
//

#include "ExN01DetectorConstruction.hh"
#include "ExN01SensitiveDetector.hh"

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
#include "G4PSCellFlux.hh"

#include "G4ParticleTable.hh"
#include "G4SDParticleFilter.hh"
//#include "G4VSDFilter.hh"

#include "DagSolid.hh"
#include "MBInterface.hpp"
#include "DagMC.hpp"

#include "DagSolidMaterial.hh"
#include "DagSolidTally.hh"

#include "../pyne/pyne.h"
//#include "pyne/particle.h"

using namespace moab;

 // dag_volumes
std::map<int, G4LogicalVolume*> dag_logical_volumes;


DagMC* dagmc = DagMC::instance(); // create dag instance
UWUW workflow_data;

ExN01DetectorConstruction::ExN01DetectorConstruction(std::string uwuw_file)
  :  world_volume_log(0)

{
  uwuw_filename = uwuw_file;

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

  workflow_data = UWUW(uwuw_filename);

  ;
}

ExN01DetectorConstruction::~ExN01DetectorConstruction()
{
}

G4VPhysicalVolume* ExN01DetectorConstruction::Construct()
{

  // load the material from the UW^2 library
  std::map<std::string,G4Material*> material_lib;
  material_lib = load_uwuw_materials(workflow_data);

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
      //      dag_logical_volumes.push_back(dag_vol_log);
      dag_logical_volumes[dag_idx]=dag_vol_log;

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
  tally_library = workflow_data.tally_library;

  // TrackLength Scorer
  //  std::vector<G4MultiFunctionalDetector*> detectors;

  // for regisistering particle tallies
  std::map<std::string,G4SDParticleFilter*> particle_filters;

  /* notes for andy
     register logical volume as sensitive only once
     then apply multiple primative sensitivities and registers
   */

  std::map<int,std::vector<std::string> > volume_part_map;
  std::vector<std::string> particles;

  // build map of vectors of volumes vs particle names
  for ( t_it = tally_library.begin() ; t_it != tally_library.end() ; ++t_it )
    {
      pyne::Tally scorer = t_it -> second; // current tally

      if( scorer.tally_type.find("Flux") != std::string::npos &&
	  scorer.entity_type.find("Volume") != std::string::npos )
	{
	  int vol_id = scorer.entity_id;
	  MBEntityHandle vol = dagmc->entity_by_id(3,vol_id); // convert id to eh
	  int vol_idx = dagmc->index_by_handle(vol); // get the volume index

	  particles = volume_part_map[vol_id];
	  particles.push_back(scorer.particle_name);
	  volume_part_map[vol_id] = particles;

	  // build the filters
	  /*
	  if( 0 < particle_filters.count(scorer.particle_name))
	    {
	      G4SDParticleFilter *filter = new G4SDParticleFilter(scorer.particle_name);
	      filter->add(pyne::particle::geant4(scorer.particle_name));
	      particle_filters[scorer.particle_name]=filter;
	    }
	  */
	}
    }

  // prints the maps
  std::map<int,std::vector<std::string> >::iterator it;
  for ( it = volume_part_map.begin() ; it != volume_part_map.end() ; ++it )
  {
    std::vector<std::string>::iterator str;
    for ( str = (it->second).begin() ; str != (it->second).end() ; ++str )
	  {
	    std::cout << (it->first) << " " << (*str) << std::endl;
    }
  }

  //  loop over the volume indices
  for ( it = volume_part_map.begin() ; it != volume_part_map.end() ; ++it )
    {
      // turn the idx into string
      std::cout << (it->first) << std::endl;
      std::stringstream int_to_string;
      int_to_string << (it->first);
      std::string idx_str = int_to_string.str();

      // get detector name
      std::string detector_name = "vol_"+idx_str+"_flux";

     // loop over the vector of particle types
     std::vector<std::string>::iterator str;
     for ( str = (it->second).begin() ; str != (it->second).end() ; ++str )
     {
        std::string particle_name = *str;
        // create particle filter
        G4SDParticleFilter *filter = new G4SDParticleFilter(particle_name);
        if(!pyne::particle::is_heavy_ion(particle_name))
           filter->add(pyne::particle::geant4(particle_name));
        else
           // add ion by getting z and a from pyne
           filter->addIon(pyne::nucname::znum(particle_name),
                          pyne::nucname::anum(particle_name));

        // create new detector
        G4VSensitiveDetector* detector = new
                   ExN01SensitiveDetector(detector_name+"/"+particle_name,
                                     "flux");
        //sets the sensitivity
        detector->SetFilter(filter);
        G4SDManager::GetSDMpointer()->AddNewDetector(detector);

        MBEntityHandle vol = dagmc->entity_by_id(3,it->first); // convert id to eh
        int vol_idx = dagmc->index_by_handle(vol); // get the volume index
        dag_logical_volumes[vol_idx]->SetSensitiveDetector(detector);
    }
  }
}
