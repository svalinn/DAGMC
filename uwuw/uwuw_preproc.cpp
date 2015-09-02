#include "uwuw_preproc.hpp"

#define DAG moab::DagMC::instance()

int main(int argc, char* argv[])
{
  std::string material_library_name = "test.h5m"; 
  // use uwuw to get materials
  UWUW *mat_lib = new UWUW();


  // load the dag file
  moab::ErrorCode rval = DAG->load_file("test.h5m");
  rval = DAG->init_OBBTree();

  std::cout << rval << std::endl;

  uwuw_preprocessor *uwuw_preproc = new uwuw_preprocessor();

  // load the materials only
  uwuw_preproc->material_library = mat_lib->load_pyne_materials(material_library_name,"/materials");
  std::cout << uwuw_preproc->material_library.size() << std::endl;
  //
  uwuw_preproc->get_dagmc_properties();

  //
  uwuw_preproc->process_materials();

  return 0;
}

uwuw_preprocessor::uwuw_preprocessor()
{
}

uwuw_preprocessor::~uwuw_preprocessor()
{
}

// get the dagmc properties
std::map<moab::EntityHandle, std::vector<std::string> > get_property_assignments( std::string property,
										 int dimension, 
										 std::string delimiters)
{
  std::map<moab::EntityHandle,std::vector<std::string> > prop_map; // to return

  std::vector<std::string> test_keywords;
  test_keywords.push_back(property);
  std::map<std::string, std::string> keyword_synonyms;

  // get initial sizes                               
  int num_entities = DAG->num_entities( dimension );

  // parse data from geometry                                                                             
  moab::ErrorCode rval = DAG->parse_properties( test_keywords, keyword_synonyms,delimiters.c_str());

  if (moab::MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  // loop over all cells                              
  for( int i = 1; i <= num_entities; ++i ) {
    //                                                
    std::vector<std::string> properties;

    // get cellid                                     
    moab::EntityHandle entity = DAG->entity_by_index( dimension, i );

    // get the group contents                         
    if( DAG->has_prop( entity, property ) )
      rval = DAG->prop_values(entity,property,properties);
    else
      properties.push_back("");

    // remove duplicates                              
    std::vector<std::string>::iterator it;
    it = std::unique(properties.begin(),properties.end());
    // resize vector to remove empty parts            
    properties.resize(std::distance(properties.begin(),it));

    // assign the map value
    prop_map[entity]=properties;
  }

  return prop_map;
}

// process the group names into unique material objects
void uwuw_preprocessor::process_materials()
{
  std::set<std::string> material_names;
  
  // vol prop iterator
  std::map<std::string,std::pair<std::string,std::string> > :: iterator it;
  // loop over each property to get the list of unique material names
  for ( it = volume_property_map.begin() ; it != volume_property_map.end() ; ++it ) {
    std::string grp_name = it->first;
    std::pair<std::string,std::string> mat_dens = it->second;
    std::string material = mat_dens.first;
    // add to the list of material names
    material_names.insert(grp_name);
    //    std::cout << material << std::endl;
    //    pyne::Material mat = material_library[material];
  }

  // 
  std::set<std::string> :: iterator s_it;
  for ( s_it = material_names.begin() ; s_it != material_names.end() ; ++s_it) {
    std::pair<std::string,std::string> mat_dens = volume_property_map[*s_it];
    std::string mat_name = mat_dens.first;
    if ( material_library.count(mat_name) == 0 ) {
      std::cout << "material " << mat_name << " was not found in the material library" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  return;
}

// get the tags on the dagmc model
void uwuw_preprocessor::get_dagmc_properties()
{
  std::map<moab::EntityHandle,std::vector<std::string> > material_assignments;
  material_assignments = get_property_assignments("mat",3,":/");
  std::map<moab::EntityHandle,std::vector<std::string> > density_assignments;
  density_assignments = get_property_assignments("rho",3,":");
  std::map<moab::EntityHandle,std::vector<std::string> > boundary_assignments;
  boundary_assignments = get_property_assignments("boundary",2,":");

  int num_cells = DAG->num_entities( 3 );

  std::vector<std::string> material_props;
  std::vector<std::string> density_props;

  pyne::Material material;

  double density;
  int material_number;

  // loop over all cells  
  for( int i = 1; i <= num_cells; ++i ) {

    density = 0.0;
    material_number = 0;

    int cellid = DAG->id_by_index( 3, i );
    moab::EntityHandle entity = DAG->entity_by_index( 3, i );

    material_props = material_assignments[entity];
    density_props = density_assignments[entity];

    if(material_props.size() == 0 ) {
      std::cout << "Volume " << cellid << " has no 'mat:' property " << std::endl;
      std::cout << "Please check your material assignments " << cellid << std::endl;
      exit(EXIT_FAILURE);
    }

    // check the that each cell's material assignment is unique
    if( material_props.size() > 1 ) {
      std::cout << "more than one material for volume with id " << cellid << std::endl;
      std::cout << cellid << " has the following material assignments" << std::endl;
      for ( int j = 0 ; j < material_props.size() ; j++ ) {
        std::cout << material_props[j] << std::endl;
      }
      std::cout << "Please check your material assignments " << cellid << std::endl;
      exit(EXIT_FAILURE);
    }
    // check the that each cell's density assignment is unique
    if(density_props.size() > 1) {
      std::cout << "More than one density specified for " << cellid <<std::endl;
      std::cout << cellid << " has the following density assignments" << std::endl;
      for ( int j = 0 ; j < density_props.size() ; j++ ) {
        std::cout << density_props[j] << std::endl;
      }
      std::cout << "Please check your density assignments " << cellid << std::endl;
      exit(EXIT_FAILURE);
    }

    // group name like mat:bob/rho:charlie
    std::string grp_name = "";
    // each unique pair of density material gets a new material object
    std::pair<std::string,std::string> mat_dens_pair;
    if (!density_props[0].empty()) {
      grp_name = "mat:"+material_props[0]+"/rho:"+density_props[0];
      mat_dens_pair.first = material_props[0];
      mat_dens_pair.second = density_props[0];
    } else {
      grp_name = "mat:"+material_props[0];
      mat_dens_pair.first = material_props[0];
      mat_dens_pair.second = "";
    }

    if(DAG->is_implicit_complement(entity)) {
      grp_name = "mat:Vaccuum";
      mat_dens_pair.first = "Vaccuum";
      mat_dens_pair.second = "";
    }

    // if the volume is the implicit compliment use cmmd line to set the 
    // material for it
      
    std::cout << "group = " << grp_name << std::endl;
    volume_property_map[grp_name] = mat_dens_pair;
  }

  /*
  int num_surfs = DAG->num_entities( 2 );

  std::vector<std::string> boundary_assignment; // boundary conditions for the current entity 
  // loop over all surfaces  
  for( int i = 1; i <= num_surfs; ++i ) {
    int surfid = DAG->id_by_index( 2, i );
    moab::EntityHandle entity = DAG->entity_by_index( 2, i );

    boundary_assignment = boundary_assignments[entity];
    if (boundary_assignment.size() != 1 )
      {
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
      lcadfile << "*";
    if (boundary_assignment[0].find("White") != std::string::npos )
      lcadfile << "+";

    lcadfile << surfid << std::endl;
  }

  // blankline 
  lcadfile << std::endl;

  // print materials
  lcadfile << "C materials from library" << std::endl;
  for(std::map<std::string,pyne::Material>::const_iterator it = material_library.begin() ;
      it != material_library.end() ; ++it )
    {
      pyne::Material new_material = (it->second);
      std::string material_card = new_material.mcnp();
      lcadfile << material_card;
    }

  // now do tallies 
  // loop over all cells 
  std::cout << "Tallies" << std::endl;
  int count = 1;
  for( std::map<std::string,pyne::Tally>::iterator it = tally_library.begin() ; it != tally_library.end() ; ++it ) {
    std::string tally_card = (it->second).mcnp(count,"mcnp5");
    lcadfile << tally_card;
    count++;
  }
  */
  return;
}
