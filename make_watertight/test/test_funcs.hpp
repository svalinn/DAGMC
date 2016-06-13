
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <algorithm>
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"

#include "gen.hpp"


moab::Interface *MBI();

/// struct to hold coordinates of skin edge, it's surface id, and a matched flag
struct coords_and_id {
  double x1;
  double y1;
  double z1;
  double x2;
  double y2;
  double z2;
  int  surf_id;
  bool matched;
  moab::EntityHandle vert1;
  moab::EntityHandle vert2;
};

// moves the vertex by dx, dy, dz
moab::ErrorCode move_vert( moab::EntityHandle vertex, double dx, double dy, double dz, bool verbose = false );
// moves the vertex by a random dx, dy, dz
moab::ErrorCode rand_vert_move( moab::EntityHandle vertex, double tol, bool verbose = false);

/// bumps the last vertex in the model by the x,y,z values given to the problem 
moab::ErrorCode single_vert_bump( moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z, bool verbose = false );

/// bumps the last vertex in the cylinder model in the R direction
moab::ErrorCode locked_pair_bump( moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z,  std::string root_name, bool verbose = false );

/// moves the last two verticies in the model in the same direction some distance less than the faceting tolerance
moab::ErrorCode locked_pair_bump_rand( moab::Range verts, double facet_tol,  std::string root_name, bool verbose = false );

/// selects a random pair of adjacent verticies and bumps them in x, y, and z 
moab::ErrorCode rand_locked_pair_bump( moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z,  std::string root_name, bool verbose = false );

/// selects a random pair of verticies from verts and moves them in random directions some distance less than the faceting tolerance
moab::ErrorCode rand_locked_pair_bump_rand( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false );

/// moves the last vertex in the model in a random direction by a distance less than the faceting tolerance
moab::ErrorCode rand_vert_bump( moab::Range verts, double facet_tol, std::string root_name, bool verbose = false );

/// moves the third to last and the last vertices in the model the same distance in x, y, and z
moab::ErrorCode adjplone_locked_pair_bump( moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z,  std::string root_name, bool verbose = false );

/// moves the third to last and the last verticies in the model in rand directions some distance less than the facet_tolerance
moab::ErrorCode adjplone_locked_pair_bump_rand( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false );

/// selects a random pair of adjacent verticies and bumps them the same distance in x, y, and z
moab::ErrorCode nonadj_locked_pair_bump( moab::Range verts, double bump_dist_x, double bump_dist_y, double bump_dist_z,  std::string root_name, bool verbose = false );

/// selects a random pair of adjacent verticies and bumps them along the theta direction a distance equal to the faceting tolerance
moab::ErrorCode nonadj_locked_pair_bump_rand( moab::Range verts, double facet_tol , std::string root_name, bool verbose = false );

/// appends "_mod" to the original file name and writes to a new .h5m 
moab::ErrorCode write_mod_file( std::string filename );

// used to clear all mesh data and reload the file as original
moab::ErrorCode reload_mesh(const char* filename,  moab::EntityHandle &meshset, bool debug = false);
