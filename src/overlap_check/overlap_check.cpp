
#include "overlap.hpp"

#include <set>
#include <memory>
#include "moab/Core.hpp"
#include "moab/GeomTopoTool.hpp"
#include "moab/GeomQueryTool.hpp"
#include "moab/ProgOptions.hpp"


using namespace moab;

int main(int argc, char* argv[]) {

  ProgOptions po("overlap_check: a tool that searches for overlaps in a DagMC geometry");

  std::string filename;

  po.addRequiredArg<std::string>("dag_file", "Path to DAGMC file to check", &filename);

  po.parseCommandLine(argc, argv);

  std::map<std::set<int>,std::array<double,3>> overlap_map;

  ErrorCode rval = check_file_for_overlaps(filename, overlap_map);
  MB_CHK_SET_ERR(rval, "Failure while checking for overlaps");

  std::cout << "Overlap locations found: " << overlap_map.size() << std::endl;

  for (const auto& entry : overlap_map) {
    std::set<int> overlap_vols = entry.first;
    std::array<double,3> loc = entry.second;

    std::cout << "Overlap Location: "
              << loc[0] << " "
              << loc[1] << " "
              << loc[2] << std::endl;
    std::cout << "Overlapping volumes: ";
    for (const auto& i : overlap_vols) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
