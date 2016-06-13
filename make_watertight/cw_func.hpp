#ifndef CWF_HPP
#define CWF_HPP

#include <iostream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/AdaptiveKDTree.hpp" // for merging verts
#include "moab/CartVect.hpp"


moab::Interface *MBI();
namespace cw_func {
/// checks the input mesh for watertightness. If check_topology=true, then the mesh will be checked topologically only, no tolerances allowed.
/// If check_topology = false, then the model will be checked for watertightness by proximity.
/// (i.e. so long as paired vertices are within tol of each other, the mesh will be considered watertight)
 moab::ErrorCode check_mesh_for_watertightness( moab::EntityHandle input_set, double tol, bool &sealed, bool test = false,  bool verbose = false , bool check_topology = false );

 int compare_by_coords(const void *a, const void *b);
 
 int compare_by_handle(const void *a, const void *b);
 
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

}

#endif
