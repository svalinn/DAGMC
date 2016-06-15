#include <iostream>
#include <string>
#include <map>

#include "uwuw.hpp"
#include "../pyne/pyne.h"
#include "DagMC.hpp"

class uwuw_preprocessor
{
 public:
  uwuw_preprocessor();
  ~uwuw_preprocessor();

 public:
  void get_dagmc_properties();
  void process_materials();

 public:
  // material_library
  std::map<std::string, pyne::Material> material_library;

 private:
  // map of all volume properties
  std::map<std::string,std::pair<std::string,std::string> > volume_property_map;

};
