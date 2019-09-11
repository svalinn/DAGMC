
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

  std::shared_ptr<Interface> MBI(new Core());

  moab::ErrorCode rval;

  rval = MBI->load_file(filename.c_str());
  MB_CHK_SET_ERR(rval, "Failed to load file: " << filename);

  std::shared_ptr<GeomTopoTool> GTT(new GeomTopoTool(MBI.get()));
  std::shared_ptr<GeomQueryTool> GQT(new GeomQueryTool(GTT.get()));

  rval = GTT->find_geomsets();
  MB_CHK_SET_ERR(rval, "Failed to detect geometry sets");

  std::cout << "Building OBB Trees..." << std::endl;

  rval = GTT->construct_obb_trees();
  MB_CHK_SET_ERR(rval, "Failed to build OBB trees");

  Range all_vols;
  rval = GTT->get_gsets_by_dimension(3, all_vols);
  MB_CHK_SET_ERR(rval, "Failed to get volumes from GTT");

  Range all_verts;
  rval = MBI->get_entities_by_type(0, MBVERTEX, all_verts);
  MB_CHK_SET_ERR(rval, "Failed to get all vertices");

  std::map<std::set<int>,std::array<double,3>> overlap_map;

  for (const auto& vert : all_verts) {
    double coords[3];
    rval = MBI->get_coords(&vert, 1, coords);
    MB_CHK_SET_ERR(rval, "Failed to get vertex coordinates for vert: " << vert);

    int num_found = 0;
    std::set<int> vols_found;
    for (const auto& vol : all_vols) {
      int result = 0;
      rval = GQT->point_in_volume(vol, coords, result);
      num_found += result;

      if (result == 1) { vols_found.insert(GTT->global_id(vol)); }
    }

    if (num_found > 1) {
        std::array<double,3> c;
        c[0] = coords[0]; c[1] = coords[1]; c[2] = coords[2];
        overlap_map[vols_found] = c;
    }
  }

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
