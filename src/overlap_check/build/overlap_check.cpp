
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
  int points_per_tri_edge = 0;

  po.addRequiredArg<std::string>("dag_file", "Path to DAGMC file to check", &filename);

  po.addOpt<int>("points per edge,p","Number of evenly-spaced points to test on each triangle edge", &points_per_tri_edge);

  po.parseCommandLine(argc, argv);

  // Load the file
  std::shared_ptr<Interface> MBI(new Core());

  ErrorCode rval = MBI->load_file(filename.c_str());
  MB_CHK_SET_ERR(rval, "Failed to load file: " << filename);

  if (points_per_tri_edge == 0) {
    std::cout << "NOTICE: " << "\n";
    std::cout << "\t Performing overlap check using triangle vertex locations only."
              << "\n"
              << "\t Use the '-p' option to check more points on the triangle edges."
              << "\n"
              << "\t Run '$dagmc_overlap_check --help' for more information."
              << "\n\n";
  }

  std::cout << "Running overlap check:" << std::endl;

  // check for overlaps
  OverlapMap overlap_map;
  rval = check_instance_for_overlaps(MBI, overlap_map, points_per_tri_edge);
  MB_CHK_SET_ERR(rval, "Failure while checking for overlaps");

  // if any overlaps are found, report them
  if (overlap_map.size() > 0) {
    report_overlaps(overlap_map);
  } else {
    std::cout << "No overlaps were found." << std::endl;
  }

  return 0;
}
