#include "dagmcmetadata.hpp"
#include <iostream>
#include <set>
#include <algorithm>

// constructor for metadata class
dagmcMetaData::dagmcMetaData(moab::DagMC* dag_ptr, bool verbosity)
{
  DAG = dag_ptr; // dagmc pointer
  verbose = verbosity;
  // these are the keywords that dagmc will understand
  // from groups if you need to process more
  // they should be added here
  metadata_keywords.push_back( "mat" );
  metadata_keywords.push_back( "rho" );
  metadata_keywords.push_back( "boundary" );
  metadata_keywords.push_back( "tally" );
  metadata_keywords.push_back( "importance" );

  // allow some synonyms
  keyword_synonyms[ "rho" ] = "density";
  keyword_synonyms[ "mat" ] = "material";
}

// destructor doesnt need to do anything
dagmcMetaData::~dagmcMetaData()
{
}

// load the property data from the dagmc instance
void dagmcMetaData::load_property_data()
{
  parse_material_data();
  parse_importance_data();
  parse_boundary_data();
  parse_tally_volume_data();
  parse_tally_surface_data();
  // finalise_counters();
}

// get the a given volume property on a given entity handle
std::string dagmcMetaData::get_volume_property(std::string property, moab::EntityHandle eh)
{
  std::string value = "";
  if (property == "material_density" ) {
    value = volume_material_property_data_eh[eh];
  } else if ( property == "material" ) {
    value = volume_material_data_eh[eh];
  } else if ( property == "density" ) {
    value = volume_density_data_eh[eh];
  } else if ( property == "importance" ) {
    value = volume_importance_data_eh[eh];
  } else if ( property == "tally" ) {
    value = tally_data_eh[eh];
  } else {
    std::cout << "Not a valid property for volumes" << std::endl;
  }
  return value;
}

// overloaded get_volume_property for indices and id's'
std::string dagmcMetaData::get_volume_property(std::string property, int vol, bool idx)
{
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
std::string dagmcMetaData::get_surface_property(std::string property, moab::EntityHandle eh)
{
  std::string value = "";
  if (property == "boundary" ) {
    value = surface_boundary_data_eh[eh];
  } else if (property == "tally") {
    value = tally_data_eh[eh];
  } else {
    std::cout << "Not a valid property for surfaces" << std::endl;
  }
  return value;
}

// overloaded get_surface_property for indices and id's'
std::string dagmcMetaData::get_surface_property(std::string property, int vol, bool idx)
{
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
void dagmcMetaData::parse_material_data()
{
  std::map<moab::EntityHandle,std::vector<std::string> > material_assignments;
  material_assignments = get_property_assignments("mat",3,":/",true);
  std::map<moab::EntityHandle,std::vector<std::string> > density_assignments;
  density_assignments = get_property_assignments("rho",3,":",true);

  int num_cells = DAG->num_entities( 3 );

  std::vector<std::string> material_props;
  std::vector<std::string> density_props;

  std::string implicit_complement_material = "";
  std::string implicit_complement_density = "";

  // loop over all cells
  for( int i = 1; i <= num_cells; ++i ) {

    int cellid = DAG->id_by_index( 3, i );
    moab::EntityHandle eh = DAG->entity_by_index( 3, i );

    material_props = material_assignments[eh];
    density_props = density_assignments[eh];

    // this is actually ok for a single volume, one that has the _comp tag at the end of it
    if( material_props.size() > 1 ) {
      // search the props for _comp
      std::size_t comp_found;
      int position;
      for ( int j = 0 ; j < material_props.size() ; j++ ) {
        comp_found = material_props[j].find("_comp");
        if ( comp_found != std::string::npos) {
          position = j;
          break;
        }
      }
      if ( comp_found != std::string::npos ) {
        // success found the _comp tag for the impl_compl material
        // set the impl_comp material for use later
        implicit_complement_material = material_props[position].substr(0,material_props[position].size()-5);
        implicit_complement_density = density_props[0];
        material_props.erase(material_props.begin()+position);
      } else {
        // failure a volume can only have a single material associated with it
        std::cout << "more than one material for volume with id " << cellid << std::endl;
        std::cout << "that does not the the _comp tag associated with it" << std::endl;
        std::cout << cellid << " has the following material assignments" << std::endl;
        for ( int j = 0 ; j < material_props.size() ; j++ ) {
          std::cout << material_props[j] << std::endl;
        }
        std::cout << "Please check your material assignments " << cellid << std::endl;
        exit(EXIT_FAILURE);
      }
    } else {
      // because of how the data are inserted into the map, ther are always
      // at least one entry, "" is nothing is found
      // if there is no material property - not failure for impl_comp
      if(material_props[0] == "" && !(DAG->is_implicit_complement(eh))) {
        std::cout << "No material property found for volume with ID " << cellid << std::endl;
        std::cout << "Every volume must have only one mat: property" << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    // this is never ok for a volume to have more than one proprety for density
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

    // check to see if the simplied naming scheme is used, by try to convert the
    // material property to an int
    if(try_to_make_int(material_props[0]) && density_props[0].empty() && !(DAG->is_implicit_complement(eh))) {
      std::cout << "Using the simplified nameing scheme without a density" << std::endl;
      std::cout << "property is forbidden, please rename the group mat:" << material_props[0] << std::endl;
      exit(EXIT_FAILURE);
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
    else if (grp_name.find("Vacuum") != std::string::npos) {
      volume_material_property_data_eh[eh] = "mat:Vacuum";
      volume_material_data_eh[eh] = "Vacuum";
    }
    // implicit complement
    else if (DAG->is_implicit_complement(eh)) {
      if(implicit_complement_material == "" ) {
        std::cout << "Implicit Complement assumed to be Vacuum" << std::endl;
        volume_material_property_data_eh[eh] = "mat:Vacuum";
        volume_material_data_eh[eh] = "Vacuum";
      } else {
        volume_material_property_data_eh[eh] = "mat:"+implicit_complement_material;
        if(implicit_complement_density != "") volume_material_property_data_eh[eh] += "/rho"+implicit_complement_density;
        volume_material_data_eh[eh] = implicit_complement_material;
        volume_density_data_eh[eh] = implicit_complement_density;
      }
    }
  }
}

// parse the importance data from the file
void dagmcMetaData::parse_importance_data()
{
  std::map<moab::EntityHandle,std::vector<std::string> > importance_assignments;
  importance_assignments = get_property_assignments("importance",3,":");

  int num_vols = DAG->num_entities( 3 );

  std::vector<std::string> importance_assignment; // importance conditions for the current entity
  // loop over all volumes set up the generic importance data
  for( int i = 1; i <= num_vols; ++i ) {
    int volid = DAG->id_by_index( 3, i );
    moab::EntityHandle eh = DAG->entity_by_index( 3, i );

    // vector of importance values
    importance_assignment = importance_assignments[eh];

    // set the value of each string for
    std::string importances = "|";
    for ( int j = 0 ; j < importance_assignment.size() ; j++ ) {
      if(importance_assignment[j] == "" ) break;
      // delimit each particle/value pair with a pipe symbol
      importances += importance_assignment[j]+"|";
      // also split to get key-value
      std::pair<std::string,std::string> pair = split_string(importance_assignment[j],"/");
      // add to the unique collection of particle names
      imp_particles.insert(pair.first);
      // insert into map too
      if(importance_map[eh].count(pair.first) == 0 ) {
        importance_map[eh][pair.first] = atof(pair.second.c_str());
      } else {
        std::cout << "Volume with ID " << volid << " has more than one importance " << std::endl;
        std::cout << "Assigned for particle type " << pair.first << std::endl;
        std::cout << "Only one importance value per volume per particle type is allowed" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    volume_importance_data_eh[eh] = importances;
  }

  // now find which regions do not have all particle importances
  // and give them importance 1.0;
  for(int i = 1 ; i <= num_vols ; ++i ) {
    std::set<std::string>::iterator it;
    for( it = imp_particles.begin() ; it != imp_particles.end() ; ++it ) {
      std::string particle = *it;
      moab::EntityHandle eh = DAG->entity_by_index( 3, i );
      if( importance_map[eh].count(particle) == 0 ) {
        if(verbose) {
          std::cout << "Warning: Volume with ID " << DAG->id_by_index(3,i);
          std::cout << " does not have an importance set for particle ";
          std::cout << particle << " assuming importance 1.0 " << std::endl;
        }
        // give this particle default importance
        importance_map[eh][particle] = 1.0;
      }
    }
  }
}

// parse the tally data from the file
void dagmcMetaData::parse_tally_volume_data()
{
  std::map<moab::EntityHandle,std::vector<std::string> > tally_assignments;
  tally_assignments = get_property_assignments("tally",3,":");

  int num_vols = DAG->num_entities( 3 );

  std::vector<std::string> tally_assignment; // tally assignments for the current entity
  // loop over all volumes
  for( int i = 1; i <= num_vols; ++i ) {
    int volid = DAG->id_by_index( 3, i );
    moab::EntityHandle eh = DAG->entity_by_index( 3, i );

    // vector of tally values
    tally_assignment = tally_assignments[eh];
    // set the value of each string for
    std::string tally = "|";
    for ( int j = 0 ; j < tally_assignment.size() ; j++ ) {
      // delimit each particle/value pair with a pipe symbol
      tally += tally_assignment[j]+"|";
    }
    tally_data_eh[eh] = tally;
  }
}

// parse the boundary data
void dagmcMetaData::parse_boundary_data()
{
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
      exit(EXIT_FAILURE);
    }
    // 2d entities have been tagged with the boundary condition property
    // ie. both surfaces and its members triangles,

    if(boundary_assignment[0].find("Reflecting") != std::string::npos )
      surface_boundary_data_eh[eh] = "Reflecting";
    if (boundary_assignment[0].find("White") != std::string::npos )
      surface_boundary_data_eh[eh] = "White";
  }
}

// parse the surface tally data from the file
void dagmcMetaData::parse_tally_surface_data()
{
  std::map<moab::EntityHandle,std::vector<std::string> > tally_assignments;
  tally_assignments = get_property_assignments("tally",2,":");

  int num_surfaces = DAG->num_entities( 2 );

  std::vector<std::string> tally_assignment; // surface tally assignments for the current entity
  // loop over all volumes
  for( int i = 1; i <= num_surfaces; ++i ) {
    int surfid = DAG->id_by_index( 2, i );
    moab::EntityHandle eh = DAG->entity_by_index( 2, i );

    // vector of tally values
    tally_assignment = tally_assignments[eh];
    // set the value of each string for
    std::string tally = "|";
    for ( int j = 0 ; j < tally_assignment.size() ; j++ ) {
      // delimit each particle/value pair with a pipe symbol
      tally += tally_assignment[j]+"|";
    }
    tally_data_eh[eh] = tally;
  }
}


// for a given property with a range of delimiters, and dimensionality
std::map<moab::EntityHandle,std::vector<std::string> > dagmcMetaData::get_property_assignments(std::string property,
    int dimension,
    std::string delimiters,
    bool remove_duplicates)
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

    // get the group contents
    if( DAG->has_prop( entity, property ) ) {
      rval = DAG->prop_values(entity,property,properties);

      // loop over the properties and check for any mention
      // of the 2nd delimiter, if so extract from 0
      // to second delimiter
      // by being here we already know he property exists
      // being found upto the first delimiter
      if(delimiters.size() > 1 ) {
        for ( int j = 0 ; j < properties.size() ; j++ ) {
          size_t npos = 0, first = npos;
          npos = properties[j].find(delimiters[1]);
          // extract from the match - which is either first
          // match or .length()
          properties[j] = properties[j].substr(0,npos);
        }
      }
    } else {
      properties.push_back("");
    }

    if (properties.size() > 1 )
      if (remove_duplicates)
        properties = remove_duplicate_properties(properties);
    // assign the map value
    prop_map[entity]=properties;
  }

  return prop_map;
}

// remove duplicate properties from the vector of properties
std::vector<std::string> dagmcMetaData::remove_duplicate_properties(std::vector<std::string> properties)
{
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
  // degenerate item(s) - should probably be fixed upstream eventually

  properties_set = set_remove_rich(properties_set);

  std::set<std::string>::iterator iter;

  std::vector<std::string> new_properties;
  // resize the array
  new_properties.resize(properties_set.size());

  // turn set back into vector
  std::copy(properties_set.begin(),properties_set.end(),new_properties.begin());

  return new_properties;
}

// from a given set remove any matches if they are found in order to keep the
// the information rich version. ie. if we find both neutron and neutron/1.0 keep
// the second one
std::set<std::string> dagmcMetaData::set_remove_rich(std::set<std::string> properties_set)
{
  std::set<std::string> new_set = properties_set;

  std::vector<std::set<std::string>::const_iterator> matches;
  std::set<std::string>::iterator set_it,toofar,it;

  std::set<std::string> to_delete;
  // loop over all elements in the set
  it = new_set.begin();
  while (it != new_set.end()) {
    // loop over the set trying to find similar names
    for ( set_it = new_set.begin(), toofar = new_set.end(); set_it != toofar; ++set_it)
      if ((*set_it).find(*it) != std::string::npos && (*set_it != *it)) {
        matches.push_back(it);
        matches.push_back(set_it);
      }

    // if there were more than 2 matches
    if( matches.size() > 1 ) {
      int smallest = 0;
      int len = 1e7;
      for ( int i = 0 ; i < matches.size() ; i++ ) {
        if( (*matches[i]).length() < len ) {
          smallest = i;
          len = (*matches[i]).length();
        }
      }
      // take the longest
      to_delete.insert(*matches[smallest]);
    }
    ++it;
  }

  // loop over the names to remove
  for ( it = to_delete.begin() ; it != to_delete.end() ; ++it ) {
    new_set.erase(new_set.find(*it));
  }

  // return the new set
  return new_set;
}

// unpack the packed string of the form delimeter<data>delimiter<data>delimiter into
// a vector of the form data[0],data[1] etc
std::vector<std::string> dagmcMetaData::unpack_string(std::string to_unpack, std::string delimiters)
{
  // loop through the string to unpack and return a vector of unpacked strings
  std::vector<size_t> locations;
  size_t npos = to_unpack.find(delimiters,0);
  // push back the first match
  locations.push_back(npos);
  while ( npos != std::string::npos) {
    npos = to_unpack.find(delimiters,npos+1);
    locations.push_back(npos);
    if(npos+1 == to_unpack.length()) break;
  }


  std::vector<std::string> unpacked_string;

  for ( int i = 0 ; i < locations.size() - 1  ; i++ ) {
    int length = locations[i+1]-1 - locations[i];
    std::string extract = to_unpack.substr(locations[i]+1,length);
    unpacked_string.push_back(extract);
  }

  // return the vector of strings
  return unpacked_string;
}

// from a string of the form key:property/key:property
// return the value of a desired key
std::string dagmcMetaData::return_property(std::string property_string, std::string property, std::string delimiter, bool chopped)
{
  std::string value =""; //value to return
  // first see if property exists
  std::size_t found_property = property_string.find(property);
  if( found_property != std::string::npos ) {
    // property found now pull out the data upto from found_property to / or to end of string
    // find the /
    std::size_t found_delimiter = property_string.find("/");
    // if we found delimiter
    if ( found_delimiter != std::string::npos) {
      // return string from found property to delimiter
      int str_length = found_delimiter - found_property;
      value = property_string.substr(found_property,str_length);
    } else {
      // return full property
      value = property_string.substr(found_property);
    }

    // if chopped, find "delimiter" and return the string after it
    if(chopped) {
      std::size_t found = value.find(delimiter);
      if(found != std::string::npos)
        value = value.substr(found+1);
      else
        value = "";
    }
  }
  return value;
}

std::pair<std::string,std::string> dagmcMetaData::split_string(std::string property_string, std::string delimiter)
{
  // first see if delimeter exists
  std::size_t found_delimiter = property_string.find(delimiter);
  std::string first = "";
  std::string second = "";
  // if we found delimiter
  if ( found_delimiter != std::string::npos) {
    // first match
    first = property_string.substr(0,found_delimiter);
    // second match
    int str_length = property_string.length() - found_delimiter;
    second = property_string.substr(found_delimiter+1,str_length);
  } else {
    std::cout << "Didnt find any delimiter" << std::endl;
    std::cout << "Returning empty strings" << std::endl;
  }
  std::pair<std::string, std::string> pair(first,second);
  return pair;
}

bool dagmcMetaData::try_to_make_int(std::string value)
{
  // try to convert the string value into an int
  char* end;

  int i = strtol(value.c_str(), &end, 10);
  if(*end == '\0')
    return true;
  else
    return false;
}
