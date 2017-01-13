//
//
// $Id$
//

#include "ExN01DetectorConstruction.hh"
#include "ExN01SensitiveDetector.hh"
#include "ExN01Analysis.hh"

#include "dagmcmetadata.hpp"

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
#include "moab/Interface.hpp"
#include "DagMC.hpp"

#include "DagSolidMaterial.hh"
#include "DagSolidTally.hh"

#include "../pyne/pyne.h"

moab::DagMC* dagmc = new moab::DagMC(); // create dag instance
dagmcMetaData* DMD;
// constructor
ExN01DetectorConstruction::ExN01DetectorConstruction(UWUW *uwuw_workflow_data)
  :  world_volume_log(0)
{
  workflow_data = uwuw_workflow_data;
}

// destructor
ExN01DetectorConstruction::~ExN01DetectorConstruction()
{
}

// the main method - takes the problem and loads
G4VPhysicalVolume* ExN01DetectorConstruction::Construct()
{
  // load the material from the UW^2 library
  std::map<std::string,G4Material*> material_lib;
  material_lib = load_uwuw_materials(workflow_data);

  G4VisAttributes * invis = new G4VisAttributes(G4VisAttributes::Invisible);

  //------------------------------------------------------ volumes
  // -- World Volume in which we place other volumes
  // we can probably do something like get the obb for the implicit complement
  // and use that as the size for the world vol
  G4double world_width  = 50000.0*cm;
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2.*world_width);

  G4Box* world_volume = new G4Box("world_volume_box",world_width,world_width,world_width);
  world_volume_log = new G4LogicalVolume(world_volume,material_lib["mat:Vacuum"],"world_vol_log",0,0,0);
  world_volume_log->SetVisAttributes(invis);
  G4PVPlacement* world_volume_phys = new G4PVPlacement(0,G4ThreeVector(),world_volume_log,
      "world_vol",0,false,0);

  // load the dagmc file - only to get the counts of volumes
  moab::ErrorCode rval = dagmc->load_file(workflow_data->full_filepath.c_str());
  if(rval != moab::MB_SUCCESS) {
    G4cout << "ERROR: Failed to load the DAGMC file " << uwuw_filename << G4endl;
    exit(1);
  }

  // build the trees
  rval = dagmc->init_OBBTree();
  if(rval != moab::MB_SUCCESS) {
    G4cout << "ERROR: Failed to build the OBB Tree" << G4endl;
    exit(1);
  }

  // attach a metadata instance
  DMD = new dagmcMetaData(dagmc);
  DMD->load_property_data();

  // get count of entities
  G4int num_of_objects = dagmc->num_entities(3);

  G4cout << "There are " << num_of_objects << " dag volumes" << G4endl;

  //Store a list of DagSolids, Logical Vols, and Physical Vols
  std::vector<DagSolid*> dag_volumes;
  std::vector<G4PVPlacement*> dag_physical_volumes;

  // load the properties from the metadata instance
  for ( int dag_idx = 1 ; dag_idx < num_of_objects ; dag_idx++ ) {
    G4String idx_str = _to_string(dag_idx);
    // get the MBEntity handle for the volume
    int dag_id = dagmc->id_by_index(3,dag_idx);
    moab::EntityHandle volume = dagmc->entity_by_id(3,dag_id);
    // get the material_name
    std::string mat_name = DMD->volume_material_property_data_eh[volume];

    // create new volume
    DagSolid* dag_vol = new DagSolid("vol_"+idx_str,dagmc,dag_idx);
    dag_volumes.push_back(dag_vol);
    // make new logical volume
    std::string material_name = mat_name;
    if( mat_name == "mat:Graveyard" || mat_name == "mat:Vacuum") {
      material_name = "mat:Vacuum"; 
    }
    
    G4LogicalVolume* dag_vol_log = new G4LogicalVolume(dag_vol,material_lib[material_name],
        "vol_"+idx_str+"_log",0,0,0);
    dag_logical_volumes[dag_idx]=dag_vol_log;
    // make a new physical placement
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

  if ( workflow_data->tally_library.size() == 0 ) {
    std::cout << "No Tallies found in the file, " << uwuw_filename << std::endl;
    std::cout << "Tallies not found, transport will happen, but no scores" << std::endl;
    return;
  }

  // create new histogram manager for tallies
  HistogramManager* HM = new HistogramManager();
  // load
  std::map<std::string,pyne::Tally>::iterator t_it;

  //  load tallies from the h5m file
  tally_library = workflow_data->tally_library;

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
  for ( t_it = tally_library.begin() ; t_it != tally_library.end() ; ++t_it ) {
    pyne::Tally scorer = t_it -> second; // current tally

    if( scorer.tally_type.find("Flux") != std::string::npos &&
        scorer.entity_type.find("Volume") != std::string::npos ) {
      int vol_id = scorer.entity_id;
      EntityHandle vol = dagmc->entity_by_id(3,vol_id); // convert id to eh
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
  for ( it = volume_part_map.begin() ; it != volume_part_map.end() ; ++it ) {
    std::vector<std::string>::iterator str;
    for ( str = (it->second).begin() ; str != (it->second).end() ; ++str ) {
      std::cout << (it->first) << " " << (*str) << std::endl;
    }
  }

  // build the filters

  // setup the data for the detectors histograms
  build_histogram();


  int sd_index = 0; // the number of sensitive detectors
  //  loop over the volume indices
  for ( it = volume_part_map.begin() ; it != volume_part_map.end() ; ++it ) {
    EntityHandle vol = dagmc->entity_by_id(3,it->first); // convert id to eh
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
    for ( str = particle_types.begin() ; str != particle_types.end() ; ++str ) {
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

  /*
  // add ebins for the given tally eventually
  std::vector<G4double> bin_bounds={1.000000E-14,2.000000E-13,2.929734E-13,4.291671E-13,6.286727E-13,
  			    9.209218E-13,1.349028E-12,1.976147E-12,2.894792E-12,4.240485E-12,
  			    6.211746E-12,9.099382E-12,1.332938E-11,1.952578E-11,2.860266E-11,
  			    4.189910E-11,6.137660E-11,8.990856E-11,1.317041E-10,1.929290E-10,
  			    2.826153E-10,4.139938E-10,5.315785E-10,6.250621E-10,6.825603E-10,
  			    8.336811E-10,8.764248E-10,1.125352E-9,1.444980E-9,1.855391E-9,
  			    2.382370E-9,3.059023E-9,3.927864E-9,5.043477E-9,6.475952E-9,
  			    8.315287E-9,1.067704E-8,1.370959E-8,1.760346E-8,2.260329E-8,
  			    2.902320E-8,3.726653E-8,4.785117E-8,6.144212E-8,7.889325E-8,
  			    1.013009E-7,1.300730E-7,1.670170E-7,2.144541E-7,2.753645E-7,
  			    3.535750E-7,4.539993E-7,5.829466E-7,7.485183E-7,9.611165E-7,
  			    1.234098E-6,1.363889E-6,1.507331E-6,1.584613E-6,1.665858E-6,
  			    1.841058E-6,2.034684E-6,2.248673E-6,2.485168E-6,2.612586E-6,
  			    2.746536E-6,2.863488E-6,3.035391E-6,3.354626E-6,3.707435E-6,
  			    4.097350E-6,4.307425E-6,4.528272E-6,5.004514E-6,5.530844E-6,
  			    6.267267E-6,7.101744E-6,8.047330E-6,9.118820E-6,1.033298E-5,
  			    1.170880E-5,1.326780E-5,1.503439E-5,1.703620E-5,1.930454E-5,
  			    2.133482E-5,2.187491E-5,2.357862E-5,2.417552E-5,2.478752E-5,
  			    2.605841E-5,2.808794E-5,3.182781E-5,3.430669E-5,3.517517E-5,
  			    3.606563E-5,4.086771E-5,4.630919E-5,5.247518E-5,5.656217E-5,
  			    5.946217E-5,6.251086E-5,6.737947E-5,7.635094E-5,8.651695E-5,
  			    9.803655E-5,1.110900E-4,1.167857E-4,1.227734E-4,1.290681E-4,
  			    1.356856E-4,1.426423E-4,1.499558E-4,1.576442E-4,1.616349E-4,
  			    1.657268E-4,1.699221E-4,1.742237E-4,1.831564E-4,1.925470E-4,
  			    2.024191E-4,2.127974E-4,2.237077E-4,2.351775E-4,2.472353E-4,
  			    2.599113E-4,2.732372E-4,2.801543E-4,2.872464E-4,2.945181E-4,
  			    3.019738E-4,3.096183E-4,3.174564E-4,3.337327E-4,3.508435E-4,
  			    3.688317E-4,3.877421E-4,4.076220E-4,4.285213E-4,4.504920E-4,
  			    4.735892E-4,4.978707E-4,5.104743E-4,5.233971E-4,5.366469E-4,
  			    5.502322E-4,5.784432E-4,6.081006E-4,6.392786E-4,6.720551E-4,
  			    7.065121E-4,7.427358E-4,7.808167E-4,8.208500E-4,8.629359E-4,
  			    9.071795E-4,9.536916E-4,9.616402E-4,9.778344E-4,1.002588E-3,
  			    1.053992E-3,1.108032E-3,1.164842E-3,1.194330E-3,1.224564E-3,
  			    1.287349E-3,1.353353E-3,1.422741E-3,1.495686E-3,1.533550E-3,
  			    1.572372E-3,1.612176E-3,1.652989E-3,1.737739E-3,1.826835E-3,
  			    1.873082E-3,1.920499E-3,1.969117E-3,2.018965E-3,2.122480E-3,
  			    2.231302E-3,2.268877E-3,2.306855E-3,2.345703E-3,2.365253E-3,
  			    2.385205E-3,2.425130E-3,2.465970E-3,2.592403E-3,2.725318E-3,
  			    2.865048E-3,3.011942E-3,3.088190E-3,3.166368E-3,3.246525E-3,
  			    3.328711E-3,3.499377E-3,3.678794E-3,3.867410E-3,4.065697E-3,
  			    4.274149E-3,4.493290E-3,4.607038E-3,4.723666E-3,4.843246E-3,
  			    4.965853E-3,5.220458E-3,5.488116E-3,5.769498E-3,5.915554E-3,
  			    5.915554E-3,6.218851E-3,6.376282E-3,6.537698E-3,6.592384E-3,
  			    6.647595E-3,6.703200E-3,6.872893E-3,7.046881E-3,7.225274E-3,
  			    7.408182E-3,7.595721E-3,7.788008E-3,7.985162E-3,8.187308E-3,
  			    8.394570E-3,8.607080E-3,8.824969E-3,9.048374E-3,9.277435E-3,
  			    9.512294E-3,9.753099E-3,1.000000E-2,1.025315E-2,1.051271E-2,
  			    1.077884E-2,1.105171E-2,1.133148E-2,1.161834E-2,1.191246E-2,
  			    1.221403E-2,1.252323E-2,1.284025E-2,1.316531E-2,1.349859E-2,
  			    1.384031E-2,1.419068E-2,1.454991E-2,1.491825E-2,1.529590E-2,
  			    1.568312E-2,1.608014E-2,1.648721E-2,1.690459E-2,1.733253E-2,
  			    1.777131E-2,1.822119E-2,1.868246E-2,1.915541E-2,1.964033E-2,
  			    2.000000E-2};
  */
  G4double bin_bounds[261] =       {1.000000E-14*GeV,2.000000E-13*GeV,2.929734E-13*GeV,4.291671E-13*GeV,6.286727E-13*GeV,
                                    9.209218E-13*GeV,1.349028E-12*GeV,1.976147E-12*GeV,2.894792E-12*GeV,4.240485E-12*GeV,
                                    6.211746E-12*GeV,9.099382E-12*GeV,1.332938E-11*GeV,1.952578E-11*GeV,2.860266E-11*GeV,
                                    4.189910E-11*GeV,6.137660E-11*GeV,8.990856E-11*GeV,1.317041E-10*GeV,1.929290E-10*GeV,
                                    2.826153E-10*GeV,4.139938E-10*GeV,5.315785E-10*GeV,6.250621E-10*GeV,6.825603E-10*GeV,
                                    8.336811E-10*GeV,8.764248E-10*GeV,1.125352E-9*GeV,1.444980E-9*GeV,1.855391E-9*GeV,
                                    2.382370E-9*GeV,3.059023E-9*GeV,3.927864E-9*GeV,5.043477E-9*GeV,6.475952E-9*GeV,
                                    8.315287E-9*GeV,1.067704E-8*GeV,1.370959E-8*GeV,1.760346E-8*GeV,2.260329E-8*GeV,
                                    2.902320E-8*GeV,3.726653E-8*GeV,4.785117E-8*GeV,6.144212E-8*GeV,7.889325E-8*GeV,
                                    1.013009E-7*GeV,1.300730E-7*GeV,1.670170E-7*GeV,2.144541E-7*GeV,2.753645E-7*GeV,
                                    3.535750E-7*GeV,4.539993E-7*GeV,5.829466E-7*GeV,7.485183E-7*GeV,9.611165E-7*GeV,
                                    1.234098E-6*GeV,1.363889E-6*GeV,1.507331E-6*GeV,1.584613E-6*GeV,1.665858E-6*GeV,
                                    1.841058E-6*GeV,2.034684E-6*GeV,2.248673E-6*GeV,2.485168E-6*GeV,2.612586E-6*GeV,
                                    2.746536E-6*GeV,2.863488E-6*GeV,3.035391E-6*GeV,3.354626E-6*GeV,3.707435E-6*GeV,
                                    4.097350E-6*GeV,4.307425E-6*GeV,4.528272E-6*GeV,5.004514E-6*GeV,5.530844E-6*GeV,
                                    6.267267E-6*GeV,7.101744E-6*GeV,8.047330E-6*GeV,9.118820E-6*GeV,1.033298E-5*GeV,
                                    1.170880E-5*GeV,1.326780E-5*GeV,1.503439E-5*GeV,1.703620E-5*GeV,1.930454E-5*GeV,
                                    2.133482E-5*GeV,2.187491E-5*GeV,2.357862E-5*GeV,2.417552E-5*GeV,2.478752E-5*GeV,
                                    2.605841E-5*GeV,2.808794E-5*GeV,3.182781E-5*GeV,3.430669E-5*GeV,3.517517E-5*GeV,
                                    3.606563E-5*GeV,4.086771E-5*GeV,4.630919E-5*GeV,5.247518E-5*GeV,5.656217E-5*GeV,
                                    5.946217E-5*GeV,6.251086E-5*GeV,6.737947E-5*GeV,7.635094E-5*GeV,8.651695E-5*GeV,
                                    9.803655E-5*GeV,1.110900E-4*GeV,1.167857E-4*GeV,1.227734E-4*GeV,1.290681E-4*GeV,
                                    1.356856E-4*GeV,1.426423E-4*GeV,1.499558E-4*GeV,1.576442E-4*GeV,1.616349E-4*GeV,
                                    1.657268E-4*GeV,1.699221E-4*GeV,1.742237E-4*GeV,1.831564E-4*GeV,1.925470E-4*GeV,
                                    2.024191E-4*GeV,2.127974E-4*GeV,2.237077E-4*GeV,2.351775E-4*GeV,2.472353E-4*GeV,
                                    2.599113E-4*GeV,2.732372E-4*GeV,2.801543E-4*GeV,2.872464E-4*GeV,2.945181E-4*GeV,
                                    3.019738E-4*GeV,3.096183E-4*GeV,3.174564E-4*GeV,3.337327E-4*GeV,3.508435E-4*GeV,
                                    3.688317E-4*GeV,3.877421E-4*GeV,4.076220E-4*GeV,4.285213E-4*GeV,4.504920E-4*GeV,
                                    4.735892E-4*GeV,4.978707E-4*GeV,5.104743E-4*GeV,5.233971E-4*GeV,5.366469E-4*GeV,
                                    5.502322E-4*GeV,5.784432E-4*GeV,6.081006E-4*GeV,6.392786E-4*GeV,6.720551E-4*GeV,
                                    7.065121E-4*GeV,7.427358E-4*GeV,7.808167E-4*GeV,8.208500E-4*GeV,8.629359E-4*GeV,
                                    9.071795E-4*GeV,9.536916E-4*GeV,9.616402E-4*GeV,9.778344E-4*GeV,1.002588E-3*GeV,
                                    1.053992E-3*GeV,1.108032E-3*GeV,1.164842E-3*GeV,1.194330E-3*GeV,1.224564E-3*GeV,
                                    1.287349E-3*GeV,1.353353E-3*GeV,1.422741E-3*GeV,1.495686E-3*GeV,1.533550E-3*GeV,
                                    1.572372E-3*GeV,1.612176E-3*GeV,1.652989E-3*GeV,1.737739E-3*GeV,1.826835E-3*GeV,
                                    1.873082E-3*GeV,1.920499E-3*GeV,1.969117E-3*GeV,2.018965E-3*GeV,2.122480E-3*GeV,
                                    2.231302E-3*GeV,2.268877E-3*GeV,2.306855E-3*GeV,2.345703E-3*GeV,2.365253E-3*GeV,
                                    2.385205E-3*GeV,2.425130E-3*GeV,2.465970E-3*GeV,2.592403E-3*GeV,2.725318E-3*GeV,
                                    2.865048E-3*GeV,3.011942E-3*GeV,3.088190E-3*GeV,3.166368E-3*GeV,3.246525E-3*GeV,
                                    3.328711E-3*GeV,3.499377E-3*GeV,3.678794E-3*GeV,3.867410E-3*GeV,4.065697E-3*GeV,
                                    4.274149E-3*GeV,4.493290E-3*GeV,4.607038E-3*GeV,4.723666E-3*GeV,4.843246E-3*GeV,
                                    4.965853E-3*GeV,5.220458E-3*GeV,5.488116E-3*GeV,5.769498E-3*GeV,5.915554E-3*GeV,
                                    6.050000E-3*GeV,6.218851E-3*GeV,6.376282E-3*GeV,6.537698E-3*GeV,6.592384E-3*GeV,
                                    6.647595E-3*GeV,6.703200E-3*GeV,6.872893E-3*GeV,7.046881E-3*GeV,7.225274E-3*GeV,
                                    7.408182E-3*GeV,7.595721E-3*GeV,7.788008E-3*GeV,7.985162E-3*GeV,8.187308E-3*GeV,
                                    8.394570E-3*GeV,8.607080E-3*GeV,8.824969E-3*GeV,9.048374E-3*GeV,9.277435E-3*GeV,
                                    9.512294E-3*GeV,9.753099E-3*GeV,1.000000E-2*GeV,1.025315E-2*GeV,1.051271E-2*GeV,
                                    1.077884E-2*GeV,1.105171E-2*GeV,1.133148E-2*GeV,1.161834E-2*GeV,1.191246E-2*GeV,
                                    1.221403E-2*GeV,1.252323E-2*GeV,1.284025E-2*GeV,1.316531E-2*GeV,1.349859E-2*GeV,
                                    1.384031E-2*GeV,1.419068E-2*GeV,1.454991E-2*GeV,1.491825E-2*GeV,1.529590E-2*GeV,
                                    1.568312E-2*GeV,1.608014E-2*GeV,1.648721E-2*GeV,1.690459E-2*GeV,1.733253E-2*GeV,
                                    1.777131E-2*GeV,1.822119E-2*GeV,1.868246E-2*GeV,1.915541E-2*GeV,1.964033E-2*GeV,
                                    2.000000E-2*GeV
                                   };

  // G4 expects a vector
  // but CI doesnt like using the C++11 std
  std::vector<G4double> bin_bounds_v;
  // hence this
  bin_bounds_v.assign(bin_bounds, bin_bounds + 261 );

  // create historgram
  //  analysisManager->CreateH1(tally_name,tally_name,10000,1e-11,100.);
  analysisManager->CreateH1(tally_name,tally_name,bin_bounds_v);
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
  if(particle_filters.count(particle_name) == 0 ) {
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

// as soon as we shift to c++11 or higher this should be removed
std::string ExN01DetectorConstruction::_to_string(int var)
{
  std::ostringstream outstr;
  outstr << var;
  std::string ret_string = outstr.str();
  return ret_string;
}
