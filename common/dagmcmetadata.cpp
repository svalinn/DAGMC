#include "DagMC.hpp"
#include "dagmcmetadata.hpp"
#include <iostream>
#include <set>
#include <algorithm> 

dagmcMetaData::dagmcMetaData(moab::DagMC *dag_ptr) {
  DAG = dag_ptr;
 
  metadata_keywords.push_back( "mat" );
  metadata_keywords.push_back( "rho" );
  metadata_keywords.push_back( "boundary" );
  metadata_keywords.push_back( "tally" );
  metadata_keywords.push_back( "importance" );

  keyword_synonyms[ "rho" ] = "density";
  keyword_synonyms[ "mat" ] = "material";
}

// destructor doesnt need to do anything
dagmcMetaData::~dagmcMetaData() {
}

// load the property data from the dagma instance
void dagmcMetaData::load_property_data(){
  parse_material_data();
  parse_importance_data();
  parse_boundary_data();
  //  parse_tally_data();
  // finalise_counters();
}

// get the a given volume property on a given entity handle
std::string dagmcMetaData::get_volume_property(std::string property, moab::EntityHandle eh) {
  std::string value = "";
  if (property == "material_density" ) {
    value = volume_material_property_data_eh[eh];
  } else if ( property == "material" ) {
    value = volume_material_data_eh[eh];
  } else if ( property == "density" ) {
    value = volume_density_data_eh[eh];
  } else if ( property == "importance" ) {
    value = volume_importance_data_eh[eh];
  } else {
    std::cout << "Not a valid property for volumes" << std::endl;
  }
  return value;
}

// overloaded get_volume_property for indices and id's'
std::string dagmcMetaData::get_volume_property(std::string property, int vol, bool idx) {
  // if this is an index query
  moab::EntityHandle eh;
  if ( idx == true) {
     eh = DAG->entity_by_index( 3, vol );
  } else {
     eh = DAG->entity_by_id( 3, vol );
  }
  return get_volume_property(property,eh);
}

// Get a given property on a surface
std::string dagmcMetaData::get_surface_property(std::string property, moab::EntityHandle eh) {
  std::string value = "";
  if (property == "material_density" ) {
    value = surface_boundary_data_eh[eh];
  } else {
    std::cout << "Not a valid property for surfaces" << std::endl;
  }
  return value;
}

// overloaded get_surface_property for indices and id's'
std::string dagmcMetaData::get_surface_property(std::string property, int vol, bool idx) {
  // if this is an index query
  moab::EntityHandle eh;
  if ( idx == true) {
     eh = DAG->entity_by_index( 2, vol );
  } else {
     eh = DAG->entity_by_id( 2, vol );
  }
  return get_surface_property(property,eh);
}


// parse the material data
void dagmcMetaData::parse_material_data() {
  std::map<moab::EntityHandle,std::vector<std::string> > material_assignments;
  material_assignments = get_property_assignments("mat",3,":/");
  std::map<moab::EntityHandle,std::vector<std::string> > density_assignments;
  density_assignments = get_property_assignments("rho",3,":");

  int num_cells = DAG->num_entities( 3 );

  std::vector<std::string> material_props;
  std::vector<std::string> density_props;

  // loop over all cells
  for( int i = 1; i <= num_cells; ++i ) {

    int cellid = DAG->id_by_index( 3, i );
    moab::EntityHandle eh = DAG->entity_by_index( 3, i );

    material_props = material_assignments[eh];
    density_props = density_assignments[eh];

    if( material_props.size() > 1 ) {
      std::cout << "more than one material for volume with id " << cellid << std::endl;
      std::cout << cellid << " has the following material assignments" << std::endl;
      for ( int j = 0 ; j < material_props.size() ; j++ ) {
        std::cout << material_props[j] << std::endl;
      }
      std::cout << "Please check your material assignments " << cellid << std::endl;
      exit(EXIT_FAILURE);
    }
    if(density_props.size() > 1) {
      std::cout << "More than one density specified for " << cellid <<std::endl;
      std::cout << cellid << " has the following density assignments" << std::endl;
      for ( int j = 0 ; j < density_props.size() ; j++ ) {
        std::cout << density_props[j] << std::endl;
      }
      std::cout << "Please check your density assignments " << cellid << std::endl;
      exit(EXIT_FAILURE);
    }

    std::string grp_name = "";

    // determine if we have a density property    
    if (!density_props[0].empty()) {
      grp_name = "mat:"+material_props[0]+"/rho:"+density_props[0];
      volume_density_data_eh[eh] = density_props[0]; 
    } else {
      grp_name = "mat:"+material_props[0];
      volume_density_data_eh[eh] = "";
    }

    // set the material value
    volume_material_property_data_eh[eh] = grp_name;
    
    // not graveyard or vacuum or implicit compliment
    if (grp_name.find("Graveyard") == std::string::npos && grp_name.find("Vacuum") == std::string::npos
        && !(DAG->is_implicit_complement(eh)) ) {
      // set the 
      volume_material_data_eh[eh] = material_props[0];
    }
    // found graveyard
    else if (grp_name.find("Graveyard") != std::string::npos) {
      volume_material_property_data_eh[eh] = "mat:Graveyard";
      volume_material_data_eh[eh] = "Graveyard";
    }
    // vacuum
    else if (grp_name.find("Vacuum") != std::string::npos || (DAG->is_implicit_complement(eh)) ) {
      volume_material_property_data_eh[eh] = "mat:Vacuum";
      volume_material_data_eh[eh] = "Vacuum";
    } 
  } 
}

// parse the importance data from the file
void dagmcMetaData::parse_importance_data() {
  std::map<moab::EntityHandle,std::vector<std::string> > importance_assignments;
  importance_assignments = get_property_assignments("importance",3,":");

  int num_vols = DAG->num_entities( 3 );

  std::vector<std::string> importance_assignment; // boundary conditions for the current entity
  // loop over all volumes
  for( int i = 1; i <= num_vols; ++i ) {
    int volid = DAG->id_by_index( 3, i );
    moab::EntityHandle eh = DAG->entity_by_index( 3, i );

    // vector of importance values
    importance_assignment = importance_assignments[eh];
    // set the value of each string for
    std::string importances = "|"; 
    for ( int j = 0 ; j < importance_assignment.size() ; j++ ) {
        // delimit each particle/value pair with a pipe symbol
        std::cout << importance_assignment[j] << std::endl;
        importances += importance_assignment[j]+"|";
     }
     volume_importance_data_eh[eh] = importances;
  }
}

// parse the boundary data
void dagmcMetaData::parse_boundary_data() {
  std::map<moab::EntityHandle,std::vector<std::string> > boundary_assignments;
  boundary_assignments = get_property_assignments("boundary",2,":");

  int num_surfs = DAG->num_entities( 2 );

  std::vector<std::string> boundary_assignment; // boundary conditions for the current entity
  // loop over all surfaces
  for( int i = 1; i <= num_surfs; ++i ) {
    int surfid = DAG->id_by_index( 2, i );
    moab::EntityHandle eh = DAG->entity_by_index( 2, i );

    boundary_assignment = boundary_assignments[eh];
    if (boundary_assignment.size() != 1 ) {
      std::cout << "More than one boundary conditions specified for " << surfid <<std::endl;
      std::cout << surfid << " has the following density assignments" << std::endl;
      for ( int j = 0 ; j < boundary_assignment.size() ; j++ ) {
        std::cout << boundary_assignment[j] << std::endl;
      }
      std::cout << "Please check your boundary condition assignments " << surfid << std::endl;

    }
    // 2d entities have been tagged with the boundary condition property
    // ie. both surfaces and its members triangles,

    if(boundary_assignment[0].find("Reflecting") != std::string::npos )
      surface_boundary_data_eh[eh] = "Reflecting";
    if (boundary_assignment[0].find("White") != std::string::npos )
      surface_boundary_data_eh[eh] = "White";
  }
}

// for a given property with a range of delimiters, and dimensionality 
std::map<moab::EntityHandle,std::vector<std::string> > dagmcMetaData::get_property_assignments(std::string property,
    int dimension, std::string delimiters )
{
  std::map<moab::EntityHandle,std::vector<std::string> > prop_map; // to return properties

   // get initial sizes
  int num_entities = DAG->num_entities( dimension );

  // parse data from geometry
  moab::ErrorCode rval = DAG->parse_properties( metadata_keywords, keyword_synonyms, delimiters.c_str());

  if (moab::MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  // loop over all entities
  for( int i = 1 ; i <= num_entities; ++i ) {
    //
    std::vector<std::string> properties;

    // get id
    moab::EntityHandle entity = DAG->entity_by_index( dimension, i );

    std::cout << "here" << std::endl;

    // get the group contents
    if( DAG->has_prop( entity, property ) )
      rval = DAG->prop_values(entity,property,properties);
    else
      properties.push_back("");

    if (properties.size() > 1 ) {
      properties = remove_duplicate_properties(properties);
      for ( int j = 0 ; j < properties.size() ; j++ ) {
	std::cout << entity << " " << j << " " << properties.size() << " " << property << " " << properties[j] << std::endl;
      }
    }
    // assign the map value
    prop_map[entity]=properties;
  }

  return prop_map;
}

std::vector<std::string> dagmcMetaData::remove_duplicate_properties(std::vector<std::string> properties){
  // remove duplicates
  std::set<std::string> properties_set; // set of properties (unique)
  
  std::vector<std::string>::iterator it;
  // loop over all properties and insert them into a set
  for ( it = properties.begin() ; it != properties.end() ; ++it) {
    properties_set.insert(*it);
  }
  // due to dagmc parse group names, the property, and value are two seperate
  // entries in array, i.e a tag like bob:charlie/bob, will return as charlie and charlie/bob
  // so we need to search each item for its more information rich partner and remove the 
  // degenerate item

  std::vector<std::set<std::string>::const_iterator> matches;
  std::set<std::string>::iterator set_it,toofar;
  // loop over all elements in the set 
  for ( set_it = properties_set.begin(), toofar = properties_set.end(); set_it != toofar; ++set_it)
    if ((*set_it).find(*properties_set.begin()) != std::string::npos) {
	matches.push_back(set_it);
    }
  // if there were more than 2 matches
  if( matches.size() > 1 ) {
    int biggest = 0;
    int len = 0;
    for ( int i = 0 ; i < matches.size() ; i++ ) {
      if( *(matches[i]).length() > len )
	biggest = i;
    }
    // take the longest 
  }
      
  std::vector<std::string> new_properties;
   // resize the array
  new_properties.resize(properties_set.size());
  
  // turn set back into vector
  std::copy(properties_set.begin(),properties_set.end(),new_properties.begin());

  return new_properties;

}
