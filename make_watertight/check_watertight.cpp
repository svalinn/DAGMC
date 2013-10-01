// ********************************************************************
// Brandon Smith
// August 2009

// This is a function to test DagMC-style mesh for watertightness. For
// now this will be a stand-alone code that uses MOAB. For volumes to 
// be watertight, the facet edges of each surface must be matched 
// one-to-one. By default this checks for topological watertightness.
// To instead use a geometrical tolerance between two vertices that are
// considered the same, pass in a tolerance.
//
// input:  h5m file name, tolerance(optional)
// output: list of unmatched facet edges and their parent surfaces, by volume.

// make CXXFLAGS=-g for debug
// make CXXFLAGS=-pg for profiling

/* Assume no MBEDGEs exist for fastest skinning
   For each volume:
     For each child surface:
       skin surface tris
       enter data into coords_and_id
     }
     match edges
   }
   Each surface is skinned twice, but the logic is simple and the memory handling 
   is easy.
*/   

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <algorithm>
#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "MBRange.hpp"
#include "MBSkinner.hpp"

#include "check_watertight_func.hpp"
#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"

MBInterface *MBI();

// struct to hold coordinates of skin edge, it's surface id, and a matched flag
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

/* qsort struct comparision function */
int compare_by_handle(const void *a, const void *b)
{
  struct coords_and_id *ia = (struct coords_and_id *)a;
  struct coords_and_id *ib = (struct coords_and_id *)b;
  if(ia->vert1 == ib->vert1) 
  {
    return (int)(ia->vert2 - ib->vert2);
  } 
  else 
  {
    return (int)(ia->vert1 - ib->vert1);
  }
  /* float comparison: returns negative if b > a 
     and positive if a > b. We multiplied result by 100.0
     to preserve decimal fraction */
} 

/* qsort struct comparision function */
// This is tricky because doubles always get rounded down to ints.
int compare_by_coords(const void *a, const void *b)
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
 
int main(int argc, char **argv) 
{

  // ******************************************************************
  // Load the h5m file and create tags.
  // ******************************************************************

  clock_t start_time;
  start_time = clock();
  // check input args
  
  if( argc < 2 || argc > 5) 
    {
    std::cout << "To check using topology of facet points:              " << std::endl;
    std::cout << "./check_watertight <filename> <verbose(true or false)>" << std::endl;
    std::cout << "To check using geometry tolerance of facet points:    " << std::endl;
    std::cout << "./check_watertight <filename> <verbose(true or false)> <tolerance>" << std::endl;
    return 1;
    }

  // load file and get tolerance from input argument
  MBErrorCode result;
  std::string filename = argv[1]; //set filename
  MBEntityHandle input_set;
  result = MBI()->create_meshset( MESHSET_SET, input_set ); //create handle to meshset
  if(MB_SUCCESS != result) 
    {
      return result;
    }

  result = MBI()->load_file( filename.c_str(), &input_set ); //load the file into the meshset
  if(MB_SUCCESS != result) 
    {
      // failed to load the file
      std::cout << "could not load file" << std::endl;
      return result;
    }

  double tol; // tolerance for two verts to be considered the same
  bool check_topology, verbose;

  if(2 == argc) // set topological check
    {
      std::cout << "topology check" << std::endl;
      check_topology = true;
      verbose = false;
    } 
  else if (3 == argc)  // set topological check with different tolerance
    {
      std::cout << "topology check" << std::endl;
      check_topology = true;
      const std::string verbose_string = argv[2];
      verbose = ( 0==verbose_string.compare("true") );
    } 
  else // otherwise do geometry check
    {
      std::cout << "geometry check";
      check_topology = false;
      tol = atof( argv[3] );
      std::cout<< " tolerance=" << tol << std::endl;
      const std::string verbose_string = argv[2];
      verbose = ( 0==verbose_string.compare("true") );
    }

  // replaced much of this code with a more modular version in check_watertight_func for testing purposes 
  std::set<int> leaky_surfs, leaky_vols;
  bool sealed, test;
  test=false;
  // is the order of the optional variables going to be a problem?
  // (i.e. we 'skipped' the variable test)
  result=check_watertight_func::check_model_for_watertightness( input_set, tol, sealed, test, verbose, check_topology);
  if(gen::error(MB_SUCCESS!=result, "could not check model for watertightness")) return result;

  clock_t end_time = clock();
  std::cout << (double) (end_time-start_time)/CLOCKS_PER_SEC << " seconds" << std::endl;
 
}


MBInterface* MBI() {
  static MBCore instance;
  return &instance;
}
