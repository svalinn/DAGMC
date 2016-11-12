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

// moab includes
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"

#include "CheckWatertight.hpp"

int main(int argc, char **argv)
{

  // ******************************************************************
  // Load the h5m file and create tags.
  // ******************************************************************

  static moab::Core instance;
  moab::Interface* mbi = &instance;

  clock_t start_time;
  start_time = clock();
  // check input args

  if( argc < 2 || argc > 5) {
    std::cout << "To check using topology of facet points:              " << std::endl;
    std::cout << "./check_watertight <filename> <verbose(true or false)>" << std::endl;
    std::cout << "To check using geometry tolerance of facet points:    " << std::endl;
    std::cout << "./check_watertight <filename> <verbose(true or false)> <tolerance>" << std::endl;
    return 1;
  }

  // load file and get tolerance from input argument
  moab::ErrorCode result;
  std::string filename = argv[1]; //set filename
  moab::EntityHandle input_set;
  result = mbi->create_meshset( moab::MESHSET_SET, input_set ); //create handle to meshset
  if(moab::MB_SUCCESS != result) {
    return result;
  }

  result = mbi->load_file( filename.c_str(), &input_set ); //load the file into the meshset
  if(moab::MB_SUCCESS != result) {
    // failed to load the file
    std::cout << "could not load file" << std::endl;
    return result;
  }

  double tol; // tolerance for two verts to be considered the same
  bool check_topology, verbose;

  if(2 == argc) { // set topological check
    std::cout << "topology check" << std::endl;
    check_topology = true;
    verbose = false;
  } else if (3 == argc) { // set topological check with different tolerance
    std::cout << "topology check" << std::endl;
    check_topology = true;
    const std::string verbose_string = argv[2];
    verbose = ( 0==verbose_string.compare("true") );
  } else { // otherwise do geometry check
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
  CheckWatertight cw = CheckWatertight(mbi);
  result= cw.check_mesh_for_watertightness( input_set, tol, sealed, test, verbose, check_topology);
  MB_CHK_SET_ERR(result, "could not check model for watertightness");

  clock_t end_time = clock();
  std::cout << (double) (end_time-start_time)/CLOCKS_PER_SEC << " seconds" << std::endl;

}


