
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

  Range all_tris;
  rval = MBI->get_entities_by_type(0, MBTRI, all_tris);

  int num_tris = all_tris.size();
  for (int tri_idx = 0; tri_idx < num_tris; tri_idx++) {
    EntityHandle tri = all_tris[tri_idx];

    Range tri_verts;
    rval = MBI->get_connectivity(&tri, 1, tri_verts);
    MB_CHK_SET_ERR(rval, "Failed to get triangle vertices");

    CartVect tri_coords[3];
    rval = MBI->get_coords(tri_verts, tri_coords[0].array());
    MB_CHK_SET_ERR(rval, "Failed to get triangle coordinates");

    // create locations for this triangle
    int num_per_edge = 100;

    std::vector<CartVect> locations;
    locations.push_back(tri_coords[0]);
    locations.push_back(tri_coords[1]);
    locations.push_back(tri_coords[2]);

    for (int i = 0; i < 3; i++) {
      CartVect edge = tri_coords[i < 2 ? i + 1 : 0] - tri_coords[i];
      CartVect vert = tri_coords[i];
      for (int j = 0; j < num_per_edge; j++) {
        double t = (double)j/(double)num_per_edge;
        locations.push_back(vert + t * edge);
      }
    }

    for (auto& loc : locations) {

      // move our point slightly off the vertex
      loc += dir * bump;

      std::set<int>vols_found;
      rval = check_location_for_overlap(GQT, all_vols, loc, dir, vols_found);
      MB_CHK_SET_ERR(rval, "Failed to for overlap at location " << loc);

      if (vols_found.size() > 1) {
        overlap_map[vols_found] = loc;
      }

      // move our point slightly off the vertex
      dir *= -1;
      loc += dir * 2.0 * bump;

      vols_found.clear();
      rval = check_location_for_overlap(GQT, all_vols, loc, dir, vols_found);
      MB_CHK_SET_ERR(rval, "Failed to for overlap at location " << loc);

      if (vols_found.size() > 1) {
        overlap_map[vols_found] = loc;
      }
    }
    prog_bar.set_value(100.0 * (double) tri_idx / (double) num_tris);
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
