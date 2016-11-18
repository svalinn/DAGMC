#include <set>
#include <iostream>

class dagmcMetaData {
  public:
   dagmcMetaData(moab::DagMC *DAGptr);
  ~dagmcMetaData();

   // load the dagmc properties into maps
   void load_property_data();

   // get a given property on a volume
   std::string get_volume_property(std::string property, moab::EntityHandle vol);
   std::string get_volume_property(std::string property, int vol, bool idx = true);

   // get a given property on a surface
   std::string get_surface_property(std::string property, moab::EntityHandle surface);
   std::string get_surface_property(std::string property, int surface, bool idx = true);

   std::vector<std::string> unpack_string(std::string to_unpack, std::string delimiters="|");
  
   std::string return_property(std::string property_string, std::string property, std::string delimiter = ":", bool chopped = true);

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

   std::map<moab::EntityHandle,std::vector<std::string> > get_property_assignments(std::string property,
										   int dimension, 
										   std::string delimiters,
										   bool remove_duplicates = true);

  std::vector<std::string> remove_duplicate_properties(std::vector<std::string> properties);

  std::set<std::string> set_remove_rich(std::set<std::string> properties_set);


  // public member variables
  public:
   // material property data map, mat:/density value
   std::map<moab::EntityHandle,std::string> volume_material_property_data_eh;

   // material data map, mat: value
   std::map<moab::EntityHandle,std::string> volume_material_data_eh;

   // density data map, rho: value
   std::map<moab::EntityHandle,std::string> volume_density_data_eh;

   // importance data map, importance: value 
   std::map<moab::EntityHandle, std::string> volume_importance_data_eh;

   // surface boundary data, rho: value
   std::map<moab::EntityHandle,std::string> surface_boundary_data_eh;

   // tally map
   std::map<moab::EntityHandle, std::string> tally_data_eh;

   // material density pairs
   // std::map<std::string, std::set<std::string> > material_density_pairs; 

  private:
  moab::DagMC *DAG; // Pointer to DAGMC instance
  
  std::vector< std::string > metadata_keywords;
  std::map< std::string, std::string > keyword_synonyms;
};
