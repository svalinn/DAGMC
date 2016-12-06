#include "uwuw_preprocessor.hpp"

// constructor
uwuw_preprocessor::uwuw_preprocessor(std::string material_library_filename, std::string dagmc_filename,
                                     std::string output_file,  bool verbosity, bool fatal_errors)
{
  // make new name concatenator class
  ncr = new name_concatenator();

  // make new DAGMC instance
  DAG = new moab::DagMC();

  // load the materials
  material_library = mat_lib.load_pyne_materials(material_library_filename,"/materials");

  // load the material objects
  // load the dag file
  moab::ErrorCode rval = DAG->load_file(dagmc_filename.c_str());

  // do the minimal DAGMC initialisation
  rval = DAG->setup_impl_compl();
  rval = DAG->setup_indices();

  // make a new dagmcmetadata class
  dmd = new dagmcMetaData(DAG);
  dmd->load_property_data();
  // set the output filename
  output_filename = output_file;

  // set the verbosity
  verbose = verbosity;
  fatal = fatal_errors;
}

// destructor
uwuw_preprocessor::~uwuw_preprocessor()
{
}

// write the new material library
void uwuw_preprocessor::write_uwuw_materials()
{
  std::map<std::string,pyne::Material> :: iterator it;

  // loop over the processed material library and write each one to the file
  for ( it = uwuw_material_library.begin() ; it != uwuw_material_library.end() ; ++it ) {
    // the current material
    pyne::Material mat = it->second;
    // write the material to the file
    if(verbose) {
      std::cout << "writing material, " << mat.metadata["name"].asString();
      std::cout << "writing material, " << mat.metadata["fluka_name"].asString();
      std::cout << " to file " << output_filename << std::endl;
    }
    // write the uwuw materials to the geometry
    mat.write_hdf5(output_filename,"/materials");
  }

  return;
}

// write the new material library
void uwuw_preprocessor::write_uwuw_tallies()
{
  std::list<pyne::Tally> :: iterator it;

  std::string tally_destination = "/tally"; // to avoid compiler warning

  // loop over the processed material library and write each one to the file
  for ( it = uwuw_tally_library.begin() ; it != uwuw_tally_library.end() ; ++it ) {
    // the current tally
    pyne::Tally tally = *it;
    if(verbose) {
      std::cout << "Writing tally " << tally.tally_name;
      std::cout << " to file " << output_filename << std::endl;
    }

    // write the material to the file
    tally.write_hdf5(output_filename,tally_destination);
  }

  return;
}


// creates a new material using the Material material as a basis
pyne::Material uwuw_preprocessor::create_new_material(pyne::Material material, std::string density)
{
  pyne::Material new_mat; // to return
  pyne::comp_map comp = material.comp;

  // make sure to bring metadata with us
  new_mat = pyne::Material(comp,1.0,material.density, 0.0, material.metadata);

  // use the name concatenator to make the fluka name
  std::string fluka_name = ncr->make_name_8bytes(material.metadata["name"].asString());

  std::string material_name;
  // make a new name
  if ( density != "" ) {
    material_name = "mat:"+material.metadata["name"].asString()+"/rho:"+density;
    new_mat.density = atof(density.c_str());
  } else {
    material_name = "mat:"+material.metadata["name"].asString();
  }

  if(verbose) {
    std::cout << "Making new material with name      : " << material_name << std::endl;
    std::cout << "                    with fluka_name: " << fluka_name << std::endl;
  }

  // tag the material with the unique name
  new_mat.metadata["name"] = material_name;
  new_mat.metadata["fluka_name"] = fluka_name;

  return new_mat;
}

// process the group names into unique material objects
void uwuw_preprocessor::process_materials()
{
  std::set<std::string> material_names;

  // vol prop iterator
  std::map<moab::EntityHandle,std::string> :: iterator it;
  std::map<moab::EntityHandle,std::string> volume_property_map = dmd->volume_material_property_data_eh;
  // loop over the volume property map

  for ( it = volume_property_map.begin() ; it != volume_property_map.end() ; ++it) {
    std::string group_name = it->second;
    if(verbose) {
      std::cout << "Making materiral group name ";
      std::cout << group_name << std::endl;
    }
    material_names.insert(group_name);
  }

  //
  std::set<std::string> :: iterator s_it;

  if (verbose) {
    std::cout << "Materials Present, :" << std::endl;
    std::map<std::string,pyne::Material> :: iterator m_it;
    for ( m_it = material_library.begin() ; m_it != material_library.end() ; ++m_it) {
      pyne::Material mat = m_it->second;
      std::cout << m_it->first << ", ";
    }
    std::cout << std::endl;
  }

  // loop over the unique material/density combinations
  for ( s_it = material_names.begin() ; s_it != material_names.end() ; ++s_it) {
    // material namesa are in the form mat:Name/rho:density
    // if we find / then we have a density token

    std::string mat_name = dmd->return_property(*s_it,"mat");
    std::string density  = dmd->return_property(*s_it,"rho");

    // find grave or vacuum
    std::size_t found_grave = mat_name.find("Graveyard");
    std::size_t found_vacuum = mat_name.find("Vacuum");

    // check for missing materials
    if ( found_grave == std::string::npos && found_vacuum == std::string::npos ) {
      if ( material_library.count(mat_name) == 0 ) {
        std::cout << "material " << mat_name << " was not found in the material library" << std::endl;
        exit(EXIT_FAILURE);
      }

      // make a new material object with the appropriate density & name
      pyne::Material new_material = create_new_material(material_library[mat_name],
                                    density);
      uwuw_material_library[*s_it]  = new_material;
    }
  }
  return;
}

// process and create all the tally objects
void uwuw_preprocessor::process_tallies()
{
  std::vector<std::string>::iterator it;

  // first volumes
  for ( int i = 1 ; i <= DAG->num_entities(3); i++ ) {
    // get the tally properties
    std::string prop = dmd->get_volume_property("tally",i,true);
    std::vector<std::string> tally_props = dmd->unpack_string(prop,"|");
    moab::EntityHandle eh = DAG->entity_by_id(3,i);
    if(tally_props.size() >= 1 && tally_props[0] != "" ) {
      for ( it = tally_props.begin() ; it != tally_props.end() ; ++it ) {
        tally_info tally_data = make_tally_groupname(*it,3,eh);
        tally_list.insert(tally_list.begin(),tally_data);
      }
    }
  }
  // now surfaces
  for ( int i = 1 ; i <= DAG->num_entities(2); i++ ) {
    // get the tally properties
    std::string prop = dmd->get_surface_property("tally",i,true);
    std::vector<std::string> tally_props = dmd->unpack_string(prop,"|");
    moab::EntityHandle eh = DAG->entity_by_id(2,i);
    if(tally_props.size() >= 1 && tally_props[0] != "" ) {
      for ( it = tally_props.begin() ; it != tally_props.end() ; ++it ) {
        tally_info tally_data = make_tally_groupname(*it,2,eh);
        tally_list.insert(tally_list.begin(),tally_data);
      }
    }
  }

  // loop over the cell tally collection and instanciate
  std::list<tally_info> :: iterator t_it;
  for ( t_it = tally_list.begin() ; t_it != tally_list.end() ; ++t_it) {
    tally_info tally_str = *t_it;
    if(verbose) {
      std::cout << tally_str.particle_name << " ";
      std::cout << tally_str.tally_name << " ";
      std::cout << tally_str.dimension << std::endl;
    }

    pyne::Tally new_tally = pyne::Tally(tally_str.tally_type,
                                        tally_str.particle_name,
                                        tally_str.entity_id,
                                        tally_str.entity_type,
                                        "", // entity name
                                        tally_str.tally_name, //tally name
                                        tally_str.entity_size,
                                        1.0); // normalisation

    uwuw_tally_library.push_back(new_tally);
  }
}

void uwuw_preprocessor::check_material_props(std::vector<std::string> material_props,
    std::vector<std::string> density_props, int cellid)
{

  if(material_props.size() == 0 ) {
    if( verbose || fatal) {
      std::cout << "Volume " << cellid << " has no 'mat:' property " << std::endl;
      std::cout << "Please check your material assignments " << cellid << std::endl;
    }
    if(fatal) exit(EXIT_FAILURE);
    no_props.push_back(cellid);
  }

  // check the that each cell's material assignment is unique
  if( material_props.size() > 1 ) {
    if( verbose || fatal ) {
      std::cout << "More than one material for volume with id " << cellid << std::endl;
      std::cout << cellid << " has the following material assignments" << std::endl;
      for ( int j = 0 ; j < material_props.size() ; j++ ) {
        std::cout << material_props[j] << std::endl;
      }
      std::cout << "Please check your material assignments " << cellid << std::endl;
    }
    if(fatal) exit(EXIT_FAILURE);
    multiple_props.push_back(cellid);
  }

  // check that there actually is an assignment
  if( material_props[0] ==  "" ) {
    if( verbose || fatal ) {
      std::cout << "Volume " << cellid << " has no 'mat:' property " << std::endl;
      std::cout << "Please check your material assignments " << cellid << std::endl;
    }
    if(fatal) exit(EXIT_FAILURE);
    blank_props.push_back(cellid);
  }

  // check the that each cell's density assignment is unique
  if(density_props.size() > 1) {
    if( verbose || fatal ) {
      std::cout << "More than one density specified for " << cellid <<std::endl;
      std::cout << cellid << " has the following density assignments" << std::endl;
      for ( int j = 0 ; j < density_props.size() ; j++ ) {
        std::cout << density_props[j] << std::endl;
      }
      std::cout << "Please check your density assignments " << cellid << std::endl;
    }
    if(fatal) exit(EXIT_FAILURE);
    multiple_densities.push_back(cellid);
  }
}

void uwuw_preprocessor::print_summary()
{
  std::cout << "+----------------------------------------------------" << std::endl;
  std::cout << "|      UWUW Summary                                  " << std::endl;
  std::cout << "+----------------------------------------------------" << std::endl;
  std::cout << "| Volumes with no property " << no_props.size() << std::endl;
  property_vector(no_props);
  std::cout << "+----------------------------------------------------" << std::endl;
  std::cout << "| Volumes with multiple properties " << multiple_props.size() << std::endl;
  property_vector(multiple_props);
  std::cout << "+----------------------------------------------------" << std::endl;
  std::cout << "| Volumes with blank properties " << blank_props.size() << std::endl;
  property_vector(blank_props);
  std::cout << "+----------------------------------------------------" << std::endl;
  std::cout << "| Volumes with multiple densities " << multiple_densities.size() << std::endl;
  property_vector(multiple_densities);
  std::cout << "+----------------------------------------------------" << std::endl;

  return;
}

// used for printing & debugging only
void uwuw_preprocessor::property_vector(std::vector<int> props)
{
  if(props.size() == 0) return;
  for ( int i = 0 ; i < props.size() ; i++ ) {
    if ( i == 0 ) {
      std::cout << "| " << props[i] << " ";
      continue;
    }
    if(i % 10 == 0 ) {
      std::cout  <<  props[i] << " " << std::endl << "| ";
    } else
      std::cout << props[i] << " ";
  }
  std::cout << std::endl;
  return;
}


tally_info uwuw_preprocessor::make_tally_groupname(std::string tally_props,
    int dimension,
    moab::EntityHandle entity)
{
  // tally props should be of the form Neutron and Flux or, Photon and Current etc

  // split tally_props into the particle and tally type
  // garenteed to be in the form, Particle/Tally
  std::string delimiter = "/";

  size_t position = 0;
  std::string particle; // the parts of the split string
  std::string tally_type;

  position = tally_props.find(delimiter);
  if( position != std::string::npos) {
    particle = tally_props.substr(0, position);
    tally_type = tally_props.substr(position+delimiter.length(),tally_props.length());
  }

  int entity_id = DAG->get_entity_id(entity);

  std::cout << tally_props << " " << dimension << " " << std::endl;
  std::cout << "particle: " << particle << std::endl;
  std::cout << tally_type << std::endl;
  check_tally_props(particle,tally_type,dimension,entity_id);

  double size = 0.0; // the area or volume of the entity
  moab::ErrorCode rval;
  std::string entity_type;
  if(dimension == 3) {
    rval = DAG->measure_volume(entity,size);
    entity_type = "Volume";
  }
  if(dimension == 2) {
    rval = DAG->measure_area(entity,size);
    entity_type = "Surface";
  }

  std::cout << entity_type << std::endl;

  // to convert the int to a string
  std::stringstream ss;
  ss << entity_id;
  std::string string_id = ss.str();

  // make the tally name
  // tally name is first two chars or particle name
  // first 2 chars of tally type
  // plus the entity id
  std::string tally_name = "";
  if(entity_id < 1000 )
    tally_name = particle.substr(0,2) + tally_type.substr(0,2)
                 + std::string(string_id);
  if(entity_id >= 1000 && entity_id < 10000 )
    tally_name = particle.substr(0,2) + tally_type.substr(0,1)
                 + std::string(string_id);
  if(entity_id >= 10000 && entity_id < 100000 )
    tally_name = particle.substr(0,1) + tally_type.substr(0,1)
                 + std::string(string_id);

  // NB We rely on the code writing out the tally name to deal with allignment

  std::transform(tally_name.begin(), tally_name.end(), tally_name.begin(), toupper);

  // keep list of particles and tally types
  tally_info tally_data;

  // set the structure data for the return
  tally_data.particle_name = particle;
  tally_data.tally_type = tally_type;
  tally_data.tally_name = tally_name;
  tally_data.dimension = dimension;
  tally_data.entity_type = entity_type;
  tally_data.entity_size = size;
  tally_data.entity_id = entity_id;

  return tally_data;
}

void uwuw_preprocessor::check_tally_props(std::string particle_name,
    std::string tally_type,
    int dimension,
    int entityid)
{

  // check to make sure the particle is allowed
  std::cout << "particle_name " << particle_name << std::endl;
  std::cout << "valid " << pyne::particle::is_valid(particle_name) << std::endl;
  if(!pyne::particle::is_valid(particle_name)) {
    if (dimension == 3) {
      std::cout << "Error: In tally on cell, ";
    } else if (dimension == 2 ) {
      std::cout << "Error: In tally on surface, ";
    }
    std::cout << particle_name << std::endl;
    std::cout << entityid << " , particle " << particle_name;
    std::cout << " not a valid particle type" << std::endl;
    exit(EXIT_FAILURE);
  }

  // check to make sure the particle is allowed
  if((tally_type != "Flux") && (tally_type != "Current")) {
    if (dimension == 3) {
      std::cout << "Error: In tally on cell, ";
    } else if (dimension == 2 ) {
      std::cout << "Error: In tally on surface, ";
    }
    std::cout << entityid << " , tally type " << tally_type;
    std::cout << " not a valid tally" << std::endl;
    exit(EXIT_FAILURE);
  }

  return;
}


// constructor
name_concatenator::name_concatenator()
{
}

// destructor
name_concatenator::~name_concatenator()
{
}

// make the fluka name from the string
std::string name_concatenator::make_name_8bytes(std::string name)
{
  // fluka name needs to be 8 chars long uppercase and unique
  std::string b8_name = name;
  std::transform(b8_name.begin(), b8_name.end(), b8_name.begin(), toupper);

  b8_name = extract_alpha_num(b8_name);

  // make the string 8 chars long
  if (b8_name.length() < 8 ) {
    b8_name = b8_name.substr(0,b8_name.length());
    b8_name.insert(b8_name.begin(), 8 - b8_name.length(), ' ');
  } else {
    b8_name = b8_name.substr(0,8);
  }

  // shift and increment if nessessary
  b8_name = shift_and_increment(b8_name);

  // insert into the used names set
  used_b8_names.insert(b8_name);

  return b8_name;
}

// function to extract only the letters [A-Z] and numbers [0-9]
std::string name_concatenator::extract_alpha_num(std::string name)
{
  // loop over the string accepting only ascii
  std::string :: iterator it;
  std::string str;

  // only accept upper case chars and numbers 1-9
  for ( it = name.begin() ; it != name.end() ; ++it ) {
    char val = *it;
    int ascii = int(val);
    if((ascii >= 65 && ascii < 91) || (ascii >= 48 && ascii < 58) )
      str += *it;
  }
  name = str;
  return name;
}

// Generates a new name but incremented by 1 relative to the last time
// the function was called.
std::string name_concatenator::shift_and_increment(std::string name)
{
  int count = 0;
  // the following statements in braces are intended to turn the 8 char long string
  // version of the name into a unique 8 character ID string, mainly for fluka
  // but has other utility in other codes
  // it should turn "Special Steel" into "SPECIALS" just by truncating and capitalising
  // however, on the second occurence of this it should turn the string into SPECIAL1,
  // then SPECIAL2, all until SPECI999
  // But, if the string is less than 8 chars anyway we should use the space first, i.e.
  // "Air, Dry N" should become "AIRDRYN ", then "AIRDRYN1", then "AIRDRY99" etc

  // if the name is unique just use it
  if(used_b8_names.count(name) == 0 )
    return name;

  // if the name isn't unique we will be appending and integer to the end
  // however we should preserve the name as much as possible, so if there are any
  // whitespace on the lhs of the string, shift left and increment there instead

  // while we haven't yet found a unique fluka_name
  while (used_b8_names.count(name) == 1 ) {
    std::size_t found = name.find_last_of(" ");
    if ( (count == 0 || count == 9 ) && found != std::string::npos)
      for ( int i = 0 ; i < 7 ; i++ ) {
        name[i] = name[i+1];
      }
    count++; // increment counter
    std::string int_as_string;
    int_to_string(count,int_as_string);
    if(count < 10) {
      name[7] = int_as_string[0];
    } else if(count >= 10 && count < 100 ) {
      name[6] = int_as_string[0];
      name[7] = int_as_string[1];
    }  else if(count >= 100 && count < 1000 ) {
      name[5] = int_as_string[0];
      name[6] = int_as_string[1];
      name[7] = int_as_string[2];
    }

    if(count == 1000) {
      std::cout << "Maximum limit of material increments reached" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  return name;
}

void name_concatenator::int_to_string(int convert, std::string &string)
{
  std::stringstream ss;
  ss << convert;
  string = ss.str();
  return;
}
