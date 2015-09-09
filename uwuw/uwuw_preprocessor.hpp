#include <iostream>
#include <string>
//#include <pair>
#include <list>
#include <map>

#include "uwuw.hpp"
#include "../pyne/pyne.h"
#include "DagMC.hpp"
#include "name_concatenator.hpp"

// convenience struct for the tally information
struct tally_info {
  std::string particle_name;
  std::string tally_name;
  std::string tally_type;
  int dimension;
  int entity_id;
  double entity_size; 
  std::string entity_name;
  std::string entity_type;
  moab::EntityHandle entity;
};

class uwuw_preprocessor {
  // public class functions
  public:
  // std string constructor
  uwuw_preprocessor(std::string material_library_filename, 
		    std::string dagmc_filename, 
 		    std::string output_file, bool verbose = false);

  // destructor
 ~uwuw_preprocessor();

  void get_dagmc_properties();
  void process_materials();
  void process_tallies();
  
  // write there new material library
  void write_uwuw_materials();

  // write there new material library
  void write_uwuw_tallies();
  
  // turn a string into the flukaname
  std::string make_fluka_name(std::string name); 

  // private class functions
  private:
    void check_material_props(std::vector<std::string> material_props,
			      std::vector<std::string> density_props,
			      int cellid);

    void check_tally_props(std::string particle_type, std::string tally_type,	 
			   int dimension, int cellid);


    void make_material_groupname(std::vector<std::string> material_props,
				 std::vector<std::string> density_props,
				 moab::EntityHandle entity,
				 std::string &grp_name,
				 std::pair<std::string,std::string > &mat_dens_pair);

    tally_info make_tally_groupname(std::string tally_props,
				    int dimension,
				    moab::EntityHandle entity);
    
  // public class members
  public:
    bool verbose;
  
  // private class members
  private:
    std::string output_filename;

    // map of all volume material properties
    std::map<std::string,std::pair<std::string,std::string> > volume_property_map;

    // map of all volume tally properties, since tallies are not unique to each volume
    // but unique to the particle and the tally type we index by that instead.
    std::list<tally_info> tally_list;

    // set of particle names for tallies
    std::set<std::string> particles;

    // set of tally types
    std::set<std::string> tallies;

    // tally_library
    std::list<pyne::Tally> uwuw_tally_library;

    // material_library
    std::map<std::string, pyne::Material> material_library;

    // material library to write out 
    std::map<std::string, pyne::Material> uwuw_material_library;

    // makes a new material
    pyne::Material create_new_material(pyne::Material material, std::string density);

    // unique naming class
    name_concatenator *ncr;

    // uwuw matlib for loading mat objs
    UWUW mat_lib;

};
