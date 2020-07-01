#ifndef SRC_DAGMC_DAGMCMETADATA_HPP_
#define SRC_DAGMC_DAGMCMETADATA_HPP_
#include <set>
#include <iostream>
#include "DagMC.hpp"

class dagmcMetaData {
 public:
  // Constructor
  dagmcMetaData(moab::DagMC* DAGptr,
                bool verbosity = false,
                bool require_density_present = true);

  // Destructor
  ~dagmcMetaData() = default;

  // load the dagmc properties into maps
  void load_property_data();

  // get a given property on a volume
  std::string get_volume_property(std::string property, moab::EntityHandle vol);
  // get a property for the specified volume, treats the vol parameter as
  // an index by default and as an ID if idx is false.
  std::string get_volume_property(std::string property, int vol, bool idx = true);

  // get a given property on a surface
  std::string get_surface_property(std::string property, moab::EntityHandle surface);
  // get a property for the specified surface, treats the surface parameter as
  // an index by default and as an ID if idx is false.
  std::string get_surface_property(std::string property, int surface, bool idx = true);

  // unpack the packed string of the form delimeter<data>delimiter<data>delimiter into
  // a vector of the form data[0],data[1] etc
  std::vector<std::string> unpack_string(std::string to_unpack, std::string delimiters = "|");

  // splits a string on the basis of the first delimiter it finds
  std::pair<std::string, std::string> split_string(std::string property_string, std::string delimiter);

  // from a string of the form key:property/key:property
  // return the value of a desired key
  std::string return_property(std::string property_string, std::string property, std::string delimiter = ":", bool chopped = true);

  // test to see if string is an int
  bool try_to_make_int(std::string value);

  // private member functions
 private:
  // parse the material data
  void parse_material_data();
  // parse the importance data
  void parse_importance_data();
  // parse the boundary data
  void parse_boundary_data();
  // parse the tally data
  void parse_tally_surface_data();
  // parse the tally data
  void parse_tally_volume_data();
  // finalise the count data
  void finalise_counters();

  // Parse property for entities with the specified dimension and delimiters.
  // Optionally remove duplicate property values if necessary.
  std::map<moab::EntityHandle, std::vector<std::string>>
  get_property_assignments(std::string property,
                           int dimension,
                           std::string delimiters,
                           bool remove_duplicates = true);

  // remove duplicate properties from the vector of properties
  std::vector<std::string> remove_duplicate_properties(std::vector<std::string> properties);

  // from a given set remove any matches if they are found in order to keep the
  // the information rich version. ie. if we find both neutron and neutron/1.0 keep
  // the second one
  std::set<std::string> set_remove_rich(std::set<std::string> properties_set);

  // public member variables
 public:
  // material property data map, mat:/density value
  // this is the full string in the form mat:<name>/density:<value>
  std::map<moab::EntityHandle, std::string> volume_material_property_data_eh;

  // material data map, mat:<name>
  std::map<moab::EntityHandle, std::string> volume_material_data_eh;

  // density data map, rho: value
  std::map<moab::EntityHandle, std::string> volume_density_data_eh;

  // importance data map, importance: value
  std::map<moab::EntityHandle, std::string> volume_importance_data_eh;

  // surface boundary data, rho: value
  std::map<moab::EntityHandle, std::string> surface_boundary_data_eh;

  // tally map
  std::map<moab::EntityHandle, std::string> tally_data_eh;

  // set to collect all particle types in the problem
  std::set<std::string> imp_particles;
  // map of importance data
  std::map<moab::EntityHandle, std::map<std::string, double>> importance_map;

  // private member variables
 private:
  moab::DagMC* DAG; // Pointer to DAGMC instance
  bool verbose; // Provide additional output while setting up and parsing properties
  bool require_density; // Require that all volumes have a specified density value
  std::vector<std::string> metadata_keywords; // Keywords supported by the metadata manager
  std::map<std::string, std::string> keyword_synonyms; // Keyword synonyms
  const std::string graveyard_str{"Graveyard"};
  const std::string vacuum_str{"Vacuum"};
};

#endif  // SRC_DAGMC_DAGMCMETADATA_HPP_
