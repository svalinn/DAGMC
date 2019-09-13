
#include "overlap.hpp"

#include "moab/GeomTopoTool.hpp"
#include "moab/GeomQueryTool.hpp"
#include "moab/ProgOptions.hpp"

using namespace moab;
ErrorCode
check_file_for_overlaps(std::shared_ptr<Interface> MBI,
                        OverlapMap& overlap_map) {

  std::shared_ptr<GeomTopoTool> GTT(new GeomTopoTool(MBI.get()));
  std::shared_ptr<GeomQueryTool> GQT(new GeomQueryTool(GTT.get()));

  ErrorCode rval = GTT->find_geomsets();
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

  std::cout << "Starting overlap check.." << std::endl;

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

  return MB_SUCCESS;
}
