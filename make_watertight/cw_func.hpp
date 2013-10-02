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
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBAdaptiveKDTree.hpp" // for merging verts
#include "MBCartVect.hpp"


MBInterface *MBI();
namespace cw_func {

 MBErrorCode check_mesh_for_watertightness( MBEntityHandle input_set, double tol, bool &sealed, bool test = false,  bool verbose = false , bool check_topology = false );

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
  MBEntityHandle vert1;
  MBEntityHandle vert2;
};

}

#endif
