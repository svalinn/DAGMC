#ifndef CWF_HPP
#define CWF_HPP

#include <iostream>
#include <iomanip>
#include <limits>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <algorithm>

// moab includes
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/CartVect.hpp"

#include "Arc.hpp"

class CheckWatertight
{

 public:
  CheckWatertight( moab::Interface* mbInterface) : mbi(mbInterface) {
    gen = new Arc(mbInterface);
  };

  ~CheckWatertight() {};
  Arc* gen;
  moab::Interface* mbi;
  moab::Interface* MBI() {
    return mbi;
  };
  /// checks the input mesh for watertightness. If check_topology=true, then the mesh will be checked topologically only, no tolerances allowed.
  /// If check_topology = false, then the model will be checked for watertightness by proximity.
  /// (i.e. so long as paired vertices are within tol of each other, the mesh will be considered watertight)
  moab::ErrorCode check_mesh_for_watertightness( moab::EntityHandle input_set, double tol, bool &sealed, bool test = false,  bool verbose = false , bool check_topology = false );

};

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

/* qsort struct comparision function */
inline int compare_by_handle(const void *a, const void *b)
{
  struct coords_and_id *ia = (struct coords_and_id *)a;
  struct coords_and_id *ib = (struct coords_and_id *)b;
  if(ia->vert1 == ib->vert1) {
    return (int)(ia->vert2 - ib->vert2);
  } else {
    return (int)(ia->vert1 - ib->vert1);
  }
  /* float comparison: returns negative if b > a
     and positive if a > b. We multiplied result by 100.0
     to preserve decimal fraction */
}

/* qsort struct comparision function */
// This is tricky because doubles always get rounded down to ints.
inline int compare_by_coords(const void *a, const void *b)
{
  struct coords_and_id *ia = (struct coords_and_id *)a;
  struct coords_and_id *ib = (struct coords_and_id *)b;
  if(ia->x1 == ib->x1) {
    if(ia->y1 == ib->y1) {
      if(ia->z1 == ib->z1) {
        if(ia->x2 == ib->x2) {
          if(ia->y2 == ib->y2) {
            if(ia->z2 == ib->z2) {
              return ia->surf_id - ib->surf_id;
            } else {
              return (ia->z2 > ib->z2) - (ia->z2 < ib->z2);
            }
          } else {
            return (ia->y2 > ib->y2) - (ia->y2 < ib->y2);
          }
        } else {
          return (ia->x2 > ib->x2) - (ia->x2 < ib->x2);
        }
      } else {
        return (ia->z1 > ib->z1) - (ia->z1 < ib->z1);
      }
    } else {
      return (ia->y1 > ib->y1) - (ia->y1 < ib->y1);;
    }
  } else {
    return (ia->x1 > ib->x1) - (ia->x1 < ib->x1);
  }
  /* float comparison: returns negative if b > a
     and positive if a > b. We multiplied result by 100.0
     to preserve decimal fraction */
}

#endif

