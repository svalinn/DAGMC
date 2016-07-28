#include <iostream>
#include <string>
#include <list>
#include <map>

#include "uwuw.hpp"
#include "../pyne/pyne.h"
#include "DagMC.hpp"
#include "name_concatenator.hpp"

/// convenience struct for the tally information
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

//===========================================================================//
/**
 * \class uwuw_preprocessor
 * \brief Defines the uwuw_preprocessor interface
 *
 * uwuw_preprocessor is a class which allows the C++ version of the uwuw preprocessor
 * program to be built. The class expects to be instantiated with 2 needed arguments
 * and 2 default arguments. The first argument defines the DAGMC file marked up with
 * appropriate instructions. The second argument defines the PyNE material library file,
 * which should contain at least the materials you specified in the DAGMC file, but can
 * contain as many as you wish. The third argument defines the output file, i.e. that which
 * you wish to output the UWUW data to. The last argument defined the verbosity of the output.
 * needed to load University of Wisconsin Unified Workflow data from the filename
 * pointed to. It contains the following private functions,
 *
 *   get_dagmc_properties(); - parses all the appropriate information from the DAGMC geometry file
 *   process_materials(); - takes all the material information found from get_dagmc_properties and turns
 *                        it into a PyNE material library
 *   process_tallies(); - takes all the tally information found from get_dagmc_properties and turns
 *                        it into a PyNE tally library
 *   write_uwuw_materials(); - writes all the material information into the output file that was specified in the
 *                           constructor. Note that the hdf5 datapath written to is "/materials"
 *   write_uwuw_tallies(); - writes all the tally information into the output file that was specified in the
 *                           constructor. Note that the hdf5 datapath written to is "/tally"
 *
 *   Once the class has been instanciated the material library, which can be read by the UWUW class, is
 *   written to the HDF5 datapath "/materials". In this library the material objects are stored along
 *   with the metadata such that each metadata name and fluka_name are unique. Similarly, the tallies
 *   are written to the HDF5 datapath "/tally". In this library the tally objects are stored along with
 *   all the tally information, such that each tally name is unique. If it is ever found that the material name
 *   or tally name isn't unique, please report that is a bug.
 **/
//===========================================================================//

class uwuw_preprocessor
{
  // public class functions
 public:
  /**
   * \brief Instanciates a new instance of the uwuw_preprocessor class, populating several important private
   * class member variables and other classes that are required.
   *
   * makes a new instance of the name_concatenator class, which takes care of making sure the material names are
   * unique. It populates the material_library data structure private to the uwuw_preprocessor class, loads the
   * DAGMC geometry into the DAG instance. Sets the class member variables output_filename (where the data will
   * be written to) and sets the verbosity of the program.
   *
   * \param[in]  material_library_filename, a string defining the full path of the material library.
   * \param[in]  dagmc_filename, a string defining the full path of the DAGMC file.
   * \param[in]  output_filename, a string defining the full path of the output file.
   * \param[in]  verbosity, a boolean value representing if vebrose (true) or non verbose(false) ouput
   *             is required.
   * \param[in]  fatal, a boolean value representing if fatal (true) then we exit on fatal errors, not fatal (false)
   *             we continue.
   */
  uwuw_preprocessor(std::string material_library_filename,
                    std::string dagmc_filename,
                    std::string output_file, bool verbose = false,
                    bool fatal = true); // constructor

  /**
   * \brief standard destructor
   */
  ~uwuw_preprocessor(); // destructor

  /**
   * \brief parses the properties in the DAGMC file and sets up class members for later use
   * , internal state relies entirely on class data, required DAGMC to have been instanciated
   * before it can be used.
   *
   * Populates the member variables volume_property_map and tally_list. The volume property map
   * stores all the correspondences between the group name i.e. that in the group like mat:Steel
   * and the material name Steel and its density should be used in the group name. The tally_list
   * stores all of the particle_name and tally_type combinations.
   */
  void get_dagmc_properties();

  /**
   * \brief Proceses the material data stored in the volume_property_map, determines the unique
   *  combinations, and makes new material object for each combination. Ultimately populating the
   * uwuw_material_library, i.e. the output material library, with all the material objects requested
   * by the input tags.
   */
  void process_materials();

  /**
   * \brief Proceses the tally data stored in the tally_list, looping over all the tallies in the list
   * and instanciating tally objects in the uwuw_tally_library.
   */
  void process_tallies();

  /**
   * \brief Writes out the PyNE Material objects, that are stored in the class member material_library variable
   * into the output file output_file as specified by the constructor.
   */
  void write_uwuw_materials();

  /**
   * \brief Writes out the PyNE Tally objects that are stored in the class member variable tally_list into
   * the output file output_file as specified by the constructor.
   */
  void write_uwuw_tallies();

  /**
   * \brief Print the summary of uwuw information
   *
   * \return void
   */
  void print_summary();

  // private class functions
 private:
  /**
   * \brief Checks that the material properties are; a) every volume has a mat: property, b) there is no more than
   *  one mat: property for each volume, c) that there is only one density assignment for the volume (either
   * implicitly with the lack of 'rho' property or by a defined 'rho' property). If any of the failure conditions
   * were met, the code will exit immediately.
   *
   * \param[in] material props, a vector of strings containing the 'mat:' properties
   * \param[in] density props, a vector of strings containing the 'rho:' properties
   * \param[in] cellid props, the integer id of the cell volume
   */
  void check_material_props(std::vector<std::string> material_props,
                            std::vector<std::string> density_props,
                            int cellid);
  /**
   * \brief Checks that for the tally properties on the current volume are; a)
   *  a valid PyNE particle type is requested, b) a valid PyNE tally type has been
   * requested.
   *
   * \param[in] particle_type, a string containing the particle type for the current tally
   * \param[in] tally_type, a string containing the tally type requested
   * \param[in] dimension, the integer representing the dimension of the entity being check
   * \param[in] entityid props, the integer id of the entity
   */
  void check_tally_props(std::string particle_type, std::string tally_type,
                         int dimension, int entityid);

  /**
   * \brief Given the vectors of strings of material & density properties and the current entity
   * makes the unique groupname by taking the first element of the material_props and density_props
   * vector and combines them to make a unqiue string in the form mat:<material_prop>/rho:<density_prop>
   * or mat:<material_prop> if no density prop has been set for the entity. Also, for convenience we
   * also return the pair of strings that makes up the material and density props.
   *
   * \param[in] material_props, vector of strings of the mat: props for this cell (by the point in time this function
   *            has been called the vector has only one element, which is the valid unqiue material name for this entity
   * \param[in] density_props, vector of strings of the rho: props for this cell (by the point in time this function has
   *            been called the vector may only have one element, containing the valid rho: property
   * \param[in] entity, the moab EntityHandle of the current volume being worked on
   * \param[out] grp_name, the unqiue group name for the entity
   * \param[out] mat_dens_pair, the mat: and rho: property pair for this entity
   */
  void make_material_groupname(std::vector<std::string> material_props,
                               std::vector<std::string> density_props,
                               moab::EntityHandle entity,
                               std::string &grp_name,
                               std::pair<std::string,std::string > &mat_dens_pair);

  /**
   * \brief Given the tally_props string the dimension make the unqiue tally name and insert the data into a
   * convenience structure tally_info.
   *
   * \param[in] tally_props, the tally property string, like Neutron/Flux
   * \param[in] dimension, the dimension of the entity being worked on
   * \param[in] entity, the MOAB entity handle of the entity being worked on
   *
   * \return the tally_info structure fully populated
   */
  tally_info make_tally_groupname(std::string tally_props,
                                  int dimension,
                                  moab::EntityHandle entity);
  /**
   * \brief Given the pyne::Material object pointed to and the new density, make a completely new
   * material object with the density pointed to, also create the unique fluka_name
   *
   * \param[in] material A PyNE material object which we would like to copy
   * \param[in] density  A string corresponding to the new density of the material
   *
   * \return a new material object with an updated density and unqiue fluka_name
   */
  pyne::Material create_new_material(pyne::Material material, std::string density);

  void property_vector(std::vector<int> props);

  // public class members
 public:
  bool verbose; ///< controls the verbosity of the class output
  bool fatal; ///< controls fatal error behavior

  // private class members
 private:
  std::string output_filename; ///< the output file which we are going to write out the data to
  std::map<std::string,std::pair<std::string,std::string> > volume_property_map; ///< map of grp name and props
  std::list<tally_info> tally_list; ///< unique list of tally information
  std::set<std::string> particles; ///< unqiue set of particles requested
  std::set<std::string> tallies; ///< set of tally types
  std::list<pyne::Tally> uwuw_tally_library; ///< unique list of tally objects
  std::map<std::string, pyne::Material> material_library; ///< material_library input by reading from library file
  std::map<std::string, pyne::Material> uwuw_material_library; ///< material library to write out to DAGMC file
  name_concatenator *ncr; ///< unique naming class pointer
  UWUW mat_lib; ///< static UWUW class for reading the material library

  std::vector<int> no_props; ///< list of cells with no properties
  std::vector<int> multiple_props; ///< list of cells with multiple properties
  std::vector<int> blank_props; ///< list of cells with blank properties
  std::vector<int> multiple_densities; ///< list of cells with multiple densities

};
