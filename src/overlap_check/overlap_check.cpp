
#include "overlap.hpp"

#include <set>
#include <memory>
#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"

using namespace moab;

int main(int argc, char* argv[]) {

  ProgOptions po("overlap_check: a tool that searches for overlaps in a DagMC geometry."
                 "This is currently a non-exhaustive search.");

  std::string filename;

  po.addRequiredArg<std::string>("dag_file", "Path to DAGMC file to check", &filename);

  po.parseCommandLine(argc, argv);

  // Load the file
  std::shared_ptr<Interface> MBI(new Core());

  ErrorCode rval = MBI->load_file(filename.c_str());
  MB_CHK_SET_ERR(rval, "Failed to load file: " << filename);

  // check for overlaps
  OverlapMap overlap_map;
  rval = check_file_for_overlaps(MBI, overlap_map);
  MB_CHK_SET_ERR(rval, "Failure while checking for overlaps");

  // if any overlaps are found, report them
  if (overlap_map.size() > 0) { report_overlaps(overlap_map); }
  else { std::cout << "No overlaps were found." << std::endl; }

  return 0;
}
