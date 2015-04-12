//
//
// $Id$
//

#include "ExN01DetectorConstruction.hh"
#include "ExN01SensitiveDetector.hh"
#include "ExN01Analysis.hh"

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

#include "G4GeometryManager.hh"
#include "G4ParticleTable.hh"
//

#include "DagSolid.hh"
#include "MBInterface.hpp"
#include "DagMC.hpp"

#include "DagSolidMaterial.hh"
#include "DagSolidTally.hh"

#include "../pyne/pyne.h"

//#include "pyne/particle.h"
using namespace moab;

DagMC* dagmc = DagMC::instance(); // create dag instance

ExN01DetectorConstruction::ExN01DetectorConstruction(UWUW *uwuw_workflow_data)
  :  world_volume_log(0)
{
  workflow_data = uwuw_workflow_data;
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


  G4double world_width  = 50000.0*cm;
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2.*world_width);

  G4Box* world_volume = new G4Box("world_volume_box",world_width,world_width,world_width);
  world_volume_log = new G4LogicalVolume(world_volume,material_lib["mat:Vacuum"],"world_vol_log",0,0,0);
  world_volume_log->SetVisAttributes(invis);
  G4PVPlacement* world_volume_phys = new G4PVPlacement(0,G4ThreeVector(),world_volume_log,
						      "world_vol",0,false,0);
<<<<<<< HEAD
=======

  G4cout << uwuw_filename.c_str() << G4endl;
  MBErrorCode rval = dagmc->load_file(uwuw_filename.c_str(),0);
  if(rval != MB_SUCCESS)
    {
      G4cout << "ERROR: Failed to load the DAGMC file " << uwuw_filename << G4endl;
      exit(1);
    }
 rval = dagmc->init_OBBTree();
 if(rval != MB_SUCCESS)
  {
    G4cout << "ERROR: Failed to build the OBB Tree" << G4endl;
    exit(1);
  }
>>>>>>> 125cb608b4fd82ec7afcbde81c4b06ebf82dab57

  MBErrorCode rval = dagmc->load_file(workflow_data->full_filepath.c_str(),0);
  if(rval != MB_SUCCESS)
    {
      G4cout << "ERROR: Failed to load the DAGMC file " << uwuw_filename << G4endl;
      exit(1);
    }
  rval = dagmc->init_OBBTree();
  if(rval != MB_SUCCESS)
    {
      G4cout << "ERROR: Failed to build the OBB Tree" << G4endl;
      exit(1);
    }
  
  G4int Num_of_survertices = dagmc->num_entities(2);
  G4int num_of_objects = dagmc->num_entities(3);

  G4cout << "There are " << num_of_objects << " dag volumes" << G4endl;

  // parse the dagmc props
  std::vector<std::string> dagsolid_keywords;
  dagsolid_keywords.push_back("mat");
  dagsolid_keywords.push_back("rho");
  std::map<std::string,std::string> synonymns;
  char const *delimeters = ":/";

  rval = dagmc->parse_properties(dagsolid_keywords,synonymns,delimeters);

<<<<<<< HEAD
=======
  rval = dagmc->parse_properties(dagsolid_keywords,synonymns,delimeters);

>>>>>>> 125cb608b4fd82ec7afcbde81c4b06ebf82dab57
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
      dag_logical_volumes[dag_idx]=dag_vol_log;

      /*
      pyne::Material mat = material_lib["mat:"+dag_material_name];
      if( dag_material_name.find("Vacuum") != std::string::npos )
	dag_vol_log->SetVisAttributes(invis);
      */

      // define physical volumes
      G4PVPlacement* dag_vol_phys = new G4PVPlacement(0,G4ThreeVector(0*cm,0*cm,0*cm),dag_vol_log,
						     "volume_"+idx_str+"_phys",
						      world_volume_log,false,0);
      dag_physical_volumes.push_back(dag_vol_phys);
    }


  return world_volume_phys;
}

// Constructs the tallies
void ExN01DetectorConstruction::ConstructSDandField()
{

  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

<<<<<<< HEAD
  if ( workflow_data->tally_library.size() == 0 )
=======
  // check for tallies first
  try
    {
      pyne::Tally tally;
      tally.from_hdf5(uwuw_filename,"/tally");
    }
  catch (const std::exception &except) // catch the exception from from_hdf5
>>>>>>> 125cb608b4fd82ec7afcbde81c4b06ebf82dab57
    {
      std::cout << "No Tallies found in the file, " << uwuw_filename << std::endl;
      std::cout << "Tallies not found, transport will happen, but no scores" << std::endl;
      return;
    }

  // create new histogram manager for tallies
  HistogramManager* HM = new HistogramManager();
  // load
  std::map<std::string,pyne::Tally>::iterator t_it;

  //  load tallies from the h5m file
<<<<<<< HEAD
  tally_library = workflow_data->tally_library;
=======
  tally_library = workflow_data.tally_library;
>>>>>>> 125cb608b4fd82ec7afcbde81c4b06ebf82dab57

  // TrackLength Scorer
  //  std::vector<G4MultiFunctionalDetector*> detectors;

  // for regisistering particle tallies

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

    // builds and keeps a store of particle filters
    BuildParticleFilter(scorer.particle_name);

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

  // build the filters

  // setup the data for the detectors histograms
  build_histogram();


  int sd_index = 0; // the number of sensitive detectors
  //  loop over the volume indices
  for ( it = volume_part_map.begin() ; it != volume_part_map.end() ; ++it )
    {
      MBEntityHandle vol = dagmc->entity_by_id(3,it->first); // convert id to eh
      int vol_idx = dagmc->index_by_handle(vol); // get the volume index

      // turn the idx into string
      std::stringstream int_to_string;
      int_to_string << vol_idx;
      std::string idx_str = int_to_string.str();

      // get detector name
      std::string detector_name = "vol_"+idx_str+"_flux";

     // increment the detector counte
     ++sd_index;

     // loop over the vector of particle types
     std::vector<std::string>::iterator str;
     std::vector<std::string> particle_types = (it->second);
     for ( str = particle_types.begin() ; str != particle_types.end() ; ++str )
     {
        // particle name
        std::string particle_name = *str;

        // create new detector
        G4cout << detector_name+'_'+particle_name << G4endl;


        int pdc = pyne::particle::id(particle_name);
        // if pdc = 0, is heavy ion, set pdc as nucid
        if ( pdc == 0 )
            pdc = pyne::nucname::id(particle_name);

        // create a histogram for each particle
        add_histogram_description(detector_name+"_"+particle_name);
        HM->add_histogram(sd_index,pdc);
      }

      // get detector volume
      double det_volume = dag_logical_volumes[vol_idx]->GetSolid()->GetCubicVolume();

      // create new detector
      G4VSensitiveDetector* detector = new ExN01SensitiveDetector(detector_name,
                                             "flux", sd_index,det_volume*cm3,
                                              HM);

      // May wish to filter particles here as well

      //sets the sensitivity
      // build the filters
      /*
      for ( str = particle_types.begin() ; str != particle_types.end() ; ++str )
        {
          // particle name
          std::string particle_name = *str;
          detector->SetFilter(particle_filters[particle_name]);
        }
       */
      // add the detector
      G4SDManager::GetSDMpointer()->AddNewDetector(detector);
      dag_logical_volumes[vol_idx]->SetSensitiveDetector(detector);
  }
  end_histogram();

  HM->print_histogram_collection();
  //exit(1);
}

/* initialise the histograms */
void ExN01DetectorConstruction::build_histogram()
{
  // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);

  // create base tuple
  analysisManager->CreateNtuple("DagGeant","TrackL");
  return;
}

/* add histogram for */
void ExN01DetectorConstruction::add_histogram_description(std::string tally_name)
{
  // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // add ebins for the given tally eventually

  // create historgram
  analysisManager->CreateH1(tally_name,tally_name,10000,1e-11,100.);
  analysisManager->CreateNtupleDColumn(tally_name);
  return;
}

/* end the histograms */
void ExN01DetectorConstruction::end_histogram()
{
  // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // finish tuple
  analysisManager->FinishNtuple();

  return;
}

// build the particle filters
void ExN01DetectorConstruction::BuildParticleFilter(std::string particle_name)
{
  // build filter if it doesnt exist
  if(particle_filters.count(particle_name) == 0 )
  {
   // create particle filter
   G4SDParticleFilter *filter = new G4SDParticleFilter(particle_name);
   if(!pyne::particle::is_heavy_ion(particle_name))
     filter->add(pyne::particle::geant4(particle_name));
   else
     // add ion by getting z and a from pyne
    filter->addIon(pyne::nucname::znum(particle_name),
                  pyne::nucname::anum(particle_name));

    particle_filters[particle_name]=filter;

   }
}
