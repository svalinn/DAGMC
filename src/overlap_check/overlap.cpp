
#include "overlap.hpp"
#include "ProgressBar.hpp"

#include "moab/GeomTopoTool.hpp"
#include "moab/GeomQueryTool.hpp"

using namespace moab;

ErrorCode check_location_for_overlap(std::shared_ptr<GeomQueryTool>& GQT,
                                     const Range& all_vols,
                                     const CartVect& loc,
                                     const CartVect& dir,
                                     std::set<int>& vols_found) {
  GeomTopoTool* GTT = GQT->gttool();
  ErrorCode rval;

  for (const auto& vol : all_vols) {
    int result = 0;
    rval = GQT->point_in_volume(vol, loc.array(), result, dir.array());
    MB_CHK_SET_ERR(rval, "Failed point in volume for Vol with id "
                   << GTT->global_id(vol)
                   << " at location " << loc);

    if (result == 1) { vols_found.insert(GTT->global_id(vol)); }
  }

  return MB_SUCCESS;
}

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

  Range all_verts;
  rval = MBI->get_entities_by_type(0, MBVERTEX, all_verts);
  MB_CHK_SET_ERR(rval, "Failed to get all vertices");

  std::cout << "Running overlap check:" << std::endl;

  ProgressBar prog_bar;

  Range all_vols;
  rval = GTT->get_gsets_by_dimension(3, all_vols);
  MB_CHK_SET_ERR(rval, "Failed to get volumes from GTT");


  CartVect dir(rand(), rand(), rand());
  dir.normalize();

  double bump = 1e-06;
  int num_verts = all_verts.size();
  for (int i = 0; i < num_verts; i++) {
    EntityHandle vert = all_verts[i];

    CartVect coords;
    rval = MBI->get_coords(&vert, 1, coords.array());
    MB_CHK_SET_ERR(rval, "Failed to get vertex coordinates for vert: " << vert);

    // move our point slightly off the vertex
    coords += dir * bump;

    std::set<int>vols_found;
    rval = check_location_for_overlap(GQT, all_vols, coords, dir, vols_found);
    MB_CHK_SET_ERR(rval, "Failed to for overlap at location " << coords);

    if (vols_found.size() > 1) {
        overlap_map[vols_found] = coords;
    }

    // move our point slightly off the vertex
    dir *= -1;
    coords += dir * 2.0 * bump;

    vols_found.clear();
    rval = check_location_for_overlap(GQT, all_vols, coords, dir, vols_found);
    MB_CHK_SET_ERR(rval, "Failed to for overlap at location " << coords);

    if (vols_found.size() > 1) {
        overlap_map[vols_found] = coords;
    }

    prog_bar.set_value(100.0 * (double) i / (double) num_verts);

  }

  return MB_SUCCESS;
}


void report_overlaps(const OverlapMap& overlap_map) {
  std::cout << "Overlap locations found: " << overlap_map.size() << std::endl;

  for (const auto& entry : overlap_map) {
    std::set<int> overlap_vols = entry.first;
    CartVect loc = entry.second;

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
}
