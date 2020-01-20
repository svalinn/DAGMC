

#include "DagMC.hpp"
#include <string>

int main() {

  // create a new DAGMC instance
  moab::DagMC* dag = new moab::DagMC();

  // get the version string
  std::string version_str;
  dag->version(&version_str);

  // write the version string to screen
  std::cout << "Found and compiled against DAGMC version: " << version_str << std::endl;

  return 0;
}
