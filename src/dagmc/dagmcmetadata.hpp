#ifndef SRC_DAGMC_DAGMCMETADATA_HPP_
#define SRC_DAGMC_DAGMCMETADATA_HPP_
#include <iostream>
#include <set>

#include "DagMC.hpp"
#include "logger.hpp"

class dagmcMetaData {
/**
 * @class dagmcMetaData
 * @brief A class to manage metadata in DAGMC models.
 *
 * The dagmcMetaData class provides functionality to parse and manage metadata
 * that defines properties associated with volumes and surfaces in a DAGMC model
 * (a .h5m file). These properties can include material assignments, density
 * values, importance values, boundary conditions, and tally assignments.
 *
 * The class provides functions to parse metadata from the model group meshsets
 * based on a set of keywords (e.g. "mat" or "rho"), retrieve metadata for
 * specific entities, and check the validity of the discovered metadata.
 */
 public:
  /**
   * @brief Constructs a new DagmcMetadata object.
   *
   * @param DAGptr A pointer to the DagMC instance associated with this metadata.
   * @param verbosity A boolean flag to control verbosity. If true, the metadata manager
   *                  will provide additional output while setting up and parsing properties.
   * @param density A boolean flag to control density requirement. If true, the metadata
   *                manager will require that all volumes have a specified density value.
   */
    dagmcMetaData(moab::DagMC* DAGptr, bool verbosity = false,
                bool require_density_present = true);

  /**
   * @brief Default destructor for the dagmcMetadata class.
   */
  ~dagmcMetaData() = default;

  /**
   * @brief Loads the DagMC properties into internal data maps.
   *
   * This function reads the properties from the associated DagMC instance and
   * stores them in internal data structures for efficient access.
   */
  void load_property_data();

  /**
   * @brief Retrieves a specified property for a given volume.
   *
   * @param property The name of the property to retrieve.
   * @param vol The handle of the volume for which to retrieve the property.
   *
   * @return Returns the value of the specified property for the given volume.
   */
  std::string get_volume_property(std::string property, moab::EntityHandle vol);

  // get a property for the specified volume, treats the vol parameter as
  // an index by default and as an ID if idx is false.

  /**
   * @brief Retrieves a specified property for a given volume.
   *
   * @param property The name of the property to retrieve.
   * @param vol The ID of the volume for which to retrieve the property.
   * @param idx Optional. If true, the vol parameter will be treated as an
   * index into the DAGMC data structures. If false, it will be treated as a
   * volume ID. Default is true.
   *
   * @return Returns the value of the specified property for the given volume.
   */
  std::string get_volume_property(std::string property, int vol,
                                  bool idx = true);

  /**
   * @brief Retrieves a specified property for a given surface.
   *
   * @param property The name of the property to retrieve.
   * @param surface The handle of the surface for which to retrieve the
   * property.
   *
   * @return Returns the value of the specified property for the given surface.
   */
  std::string get_surface_property(std::string property,
                                   moab::EntityHandle surface);

  /**
   * @brief Retrieves a specified property for a given surface.
   *
   * @param property The name of the property to retrieve.
   * @param surface The ID of the surface for which to retrieve the property.
   * @param idx Optional. If true, the surface parameter will be treated as an
   * index into the DAGMC data structures. If false, it will be treated as a
   * surface ID. Default is true.
   *
   * @return Returns the value of the specified property for the given surface.
   */
  std::string get_surface_property(std::string property, int surface,
                                   bool idx = true);

  /**
   * @brief Unpacks a string into a vector of substrings.
   *
   * The input string is expected to be of the form
   * "delimiter<data>delimiter<data>delimiter". The function will unpack this
   * into a vector of the form "data[0],data[1],...".
   *
   * @param to_unpack The string to unpack.
   * @param delimiters The delimiters used to separate the data in the string.
   *                   Default is "|".
   *
   * @return Returns a vector of substrings extracted from the input string.
   */
  std::vector<std::string> unpack_string(std::string to_unpack,
                                         std::string delimiters = "|");

  /**
   * @brief Splits a string on the basis of the first delimiter it finds.
   *
   * @param property_string The string to be split.
   * @param delimiter The delimiter on which to split the string.
   *
   * @return Returns a pair of strings. The first element of the pair is the
   *         substring before the first occurrence of the delimiter. The second
   *         element is the substring after the first occurrence of the
   *         delimiter.
   */
  std::pair<std::string, std::string> split_string(std::string property_string,
                                                   std::string delimiter);

  /**
   * @brief Extracts the value of a desired key from a property string.
   *
   * The property string is expected to be of the form
   * "key:property/key:property".
   *
   * @param property_string The string containing the properties.
   * @param property The key for which to return the value.
   * @param delimiter The delimiter used to separate the keys and values in the
   *                  property string. Default is ":".
   * @param chopped Optional. If true, the function will return the property
   *                value after removing the key and delimiter. If false, it
   *                will return the property value as is. Default is true.
   *
   * @return Returns the value of the specified key from the property string.
   */
  std::string return_property(std::string property_string, std::string property,
                              std::string delimiter = ":", bool chopped = true);

  /**
   * @brief Tests if a string can be converted to an integer.
   *
   * @param value The string to be tested.
   *
   * @return Returns true if the string can be converted to an integer,
   *         false otherwise.
   */
  bool try_to_make_int(std::string value);

  // private member functions
 private:
  /**
   * @brief Parses the material and density assignmens from group meshsets in
   * the associated DagMC instance using the "mat" and "rho" keywords,
   * respectively. Populates the volume_material_property_data_eh,
   * volume_material_data_eh, and volume_density_data_eh maps.
   */
  void parse_material_data();

  /**
   * @brief Parses the volume importance data from group meshsets in the
   * associated DagMC instance using the "importance" keyword. Populates the
   * volume_importance_data_eh map.
   */
  void parse_importance_data();

  /**
   * @brief Parses the boundary data from group meshsets in the the associated
   * DagMC instance using the "boundary" keyword. Populates the
   * surface_boundary_data_eh map.
   */
  void parse_boundary_data();

  /**
   * @brief Parses the tally surface data from group meshsets in the associated
   * DagMC instance using the "tally" keyword. Populates the tally_data_eh map.
   */
  void parse_tally_surface_data();

  /**
   * @brief Parses the tally volume data from group meshsets in the associated
   * DagMC instance using the "tally" keyword. Populates the tally_data_eh map.
   */
  void parse_tally_volume_data();

  /**
   * @brief Parses property for entities with the specified dimension and
   *        delimiters. Optionally removes duplicate property values if
   *        necessary.
   *
   * @param property The name of the property to retrieve.
   * @param dimension The dimension of the entities for which to retrieve the
   *                  property.
   * @param delimiters The delimiters used to separate the data in the string.
   * @param remove_duplicates Optional. If true, the function will remove
   *                          duplicate property values. Default is true.
   *
   * @return Returns a map where the keys are the entity handles and the values
   *         are vectors of property values.
   */
  std::map<moab::EntityHandle, std::vector<std::string>>
  get_property_assignments(std::string property, int dimension,
                           std::string delimiters,
                           bool remove_duplicates = true);

  /**
   * @brief Removes duplicate properties from a vector of properties.
   *
   * Group names are parsed in such a way that the property and value are two separate
   * entries in the array. For example, a tag like "particle:neutron/1.0" will return as
   * "neutron" and "neutron/1.0". This function searches each item for its more
   * information-rich partner and removes the degenerate item(s).
   *
   * @param properties The vector of properties from which to remove duplicates.
   *
   * @return Returns a vector of properties with all duplicates removed.
   */
  std::vector<std::string> remove_duplicate_properties(
      std::vector<std::string> properties);

  /**
   * @brief Removes less informative properties from a set of properties.
   *
   * If a property and a more informative version of that property (e.g.,
   * "neutron" and "neutron/1.0") are both present in the set, this function
   * removes the less informative version.
   *
   * @param properties_set The set of properties from which to remove less
   *                       informative properties.
   *
   * @return Returns a set of properties with all less informative properties
   *         removed.
   */
  std::set<std::string> set_remove_rich(std::set<std::string> properties_set);

  std::string to_lower(const std::string input);

  // public member variables
 public:
  /**
   * @brief Map storing the material property data for each volume.
   *
   * The keys are the entity handles of the volumes. The values are strings
   * in the form "mat:<name>/density:<value>".
   */
  std::map<moab::EntityHandle, std::string> volume_material_property_data_eh;

  /**
   * @brief Map storing the material data for each volume.
   *
   * The keys are the entity handles of the volumes. The values are strings
   * in the form "mat:<name>".
   */
  std::map<moab::EntityHandle, std::string> volume_material_data_eh;

  /**
   * @brief Map storing the density data for each volume.
   *
   * The keys are the entity handles of the volumes. The values are strings
   * in the form "rho:<value>".
   */
  std::map<moab::EntityHandle, std::string> volume_density_data_eh;

  /**
   * @brief Map storing the importance data for each volume.
   *
   * The keys are the entity handles of the volumes. The values are strings
   * in the form "importance:<value>".
   */
  std::map<moab::EntityHandle, std::string> volume_importance_data_eh;

  /**
   * @brief Map storing the boundary data for each surface.
   *
   * The keys are the entity handles of the surfaces. The values are strings
   * in the form "boundary:<value>".
   */
  std::map<moab::EntityHandle, std::string> surface_boundary_data_eh;

  /**
   * @brief Map storing the tally data.
   *
   * The keys are the entity handles. The values are strings representing
   * the tally data.
   */
  std::map<moab::EntityHandle, std::string> tally_data_eh;

  /**
   * @brief Set to collect all particle types in the problem.
   */
  std::set<std::string> imp_particles;

  /**
   * @brief Map storing the importance data for each entity.
   *
   * The keys are the entity handles. The values are maps where the keys are
   * particle types and the values are the importance values for those
   * particles.
   */
  std::map<moab::EntityHandle, std::map<std::string, double>> importance_map;

  // private member variables
 private:
  /**
   * @brief Pointer to the DAGMC instance.
   */
  moab::DagMC* DAG;

  /**
   * @brief Flag to control verbosity.
   *
   * If true, the metadata manager will provide additional output while setting
   * up and parsing properties.
   */
  bool verbose;

  /**
   * @brief Flag to control density requirement.
   *
   * If true, the metadata manager will require that all volumes have a
   * specified density value.
   */
  bool require_density;

  /**
   * @brief List of keywords supported by the metadata manager.
   */
  std::vector<std::string> metadata_keywords;

  /**
   * @brief Map of property keyword synonyms.
   *
   */
  std::map<std::string, std::string> keyword_synonyms;

  // Some constant keyword values
  const std::string graveyard_str{"Graveyard"};
  const std::string vacuum_str{"Vacuum"};
  const std::string reflecting_str{"Reflecting"};
  const std::string white_str{"White"};
  const std::string periodic_str{"Periodic"};

  DagMC_Logger logger;
};

#endif  // SRC_DAGMC_DAGMCMETADATA_HPP_
