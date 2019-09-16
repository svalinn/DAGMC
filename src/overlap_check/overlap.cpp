
#include "overlap.hpp"
#include "ProgressBar.hpp"

#include "moab/GeomTopoTool.hpp"
#include "moab/GeomQueryTool.hpp"

using namespace moab;

ErrorCode check_location_for_overlap(std::shared_ptr<GeomQueryTool>& GQT,
                                     const Range& all_vols,
                                     CartVect& loc,
                                     CartVect& dir,
                                     OverlapMap& overlap_map) {
  ErrorCode rval;

  GeomTopoTool* GTT = GQT->gttool();
  std::set<int> vols_found;
  double bump = 1E-06;

  // move the point slightly off the vertex
  loc += dir * bump;

  for (const auto& vol : all_vols) {
    int result = 0;
    rval = GQT->point_in_volume(vol, loc.array(), result, dir.array());
    MB_CHK_SET_ERR(rval, "Failed point in volume for Vol with id "
                   << GTT->global_id(vol)
                   << " at location " << loc);

    if (result == 1) {
      vols_found.insert(GTT->global_id(vol));
    }
  }

  if (vols_found.size() > 1) {
    overlap_map[vols_found] = loc;
  }

  // move the point slightly off the vertex
  dir *= -1;
  loc += dir * 2.0 * bump;
  vols_found.clear();

  for (const auto& vol : all_vols) {
    int result = 0;
    rval = GQT->point_in_volume(vol, loc.array(), result, dir.array());
    MB_CHK_SET_ERR(rval, "Failed point in volume for Vol with id "
                   << GTT->global_id(vol)
                   << " at location " << loc);

    if (result == 1) {
      vols_found.insert(GTT->global_id(vol));
    }
  }

  if (vols_found.size() > 1) {
    overlap_map[vols_found] = loc;
  }

  return MB_SUCCESS;
}

ErrorCode
check_instance_for_overlaps(std::shared_ptr<Interface> MBI,
                            OverlapMap& overlap_map,
                            int pnts_per_edge) {

  std::shared_ptr<GeomTopoTool> GTT(new GeomTopoTool(MBI.get()));
  std::shared_ptr<GeomQueryTool> GQT(new GeomQueryTool(GTT.get()));

  ErrorCode rval = GTT->find_geomsets();
  MB_CHK_SET_ERR(rval, "Failed to detect geometry sets");

  rval = GTT->construct_obb_trees();
  MB_CHK_SET_ERR(rval, "Failed to build OBB trees");

  Range all_verts;
  rval = MBI->get_entities_by_type(0, MBVERTEX, all_verts);
  MB_CHK_SET_ERR(rval, "Failed to get all vertices");

  Range all_tris;
  rval = MBI->get_entities_by_type(0, MBTRI, all_tris);
  MB_CHK_SET_ERR(rval, "Failed to get all triangles");

  Range all_edges;
  rval = MBI->get_adjacencies(all_tris, 1, true, all_edges, Interface::UNION);
  MB_CHK_SET_ERR(rval, "Failed to get triangle edges");

  Range all_vols;
  rval = GTT->get_gsets_by_dimension(3, all_vols);
  MB_CHK_SET_ERR(rval, "Failed to get volumes from GTT");

  // number of locations we'll be checking
  int num_locations = all_verts.size() + pnts_per_edge * all_edges.size();
  int num_checked = 0;

  CartVect dir(rand(), rand(), rand());
  dir.normalize();

  ProgressBar prog_bar;

  // first check all triangle vertex locations
  for (const auto& vert : all_verts) {
    CartVect loc;
    rval = MBI->get_coords(&vert, 1, loc.array());

    rval = check_location_for_overlap(GQT, all_vols, loc, dir, overlap_map);
    MB_CHK_SET_ERR(rval, "Failed to check point for overlap");

    prog_bar.set_value(100.0 * (double) num_checked++ / (double) num_locations);
  }

  // if we aren't checking along edges, return
  if (pnts_per_edge == 0) {
    return MB_SUCCESS;
  }

  // now check along triangle edges
  // (curve edges are likely in here too,
  //  but it isn't hurting anything to check more locations)
  for (const auto& edge : all_edges) {
    Range edge_verts;
    rval = MBI->get_connectivity(&edge, 1, edge_verts);
    MB_CHK_SET_ERR(rval, "Failed to get triangle vertices");

    CartVect edge_coords[2];
    rval = MBI->get_coords(edge_verts, edge_coords[0].array());
    MB_CHK_SET_ERR(rval, "Failed to get triangle coordinates");

    std::vector<CartVect> locations;

    CartVect edge_vec = edge_coords[1] - edge_coords[0];
    CartVect& start = edge_coords[0];

    // create locations along the edge to check
    for (int j = 1; j <= pnts_per_edge; j++) {
      double t = (double)j / (double)pnts_per_edge;
      locations.push_back(start + t * edge_vec);
    }

    // check edge locations
    for (auto& loc : locations) {
      rval = check_location_for_overlap(GQT, all_vols, loc, dir, overlap_map);
      MB_CHK_SET_ERR(rval, "Failed to check point for overlap");

      prog_bar.set_value(100.0 * (double) num_checked++ / (double) num_locations);
    }
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
