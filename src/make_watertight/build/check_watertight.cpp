// ********************************************************************
// Brandon Smith, August 2009
// Jonathan Shimwell, September 2019
// Patrick Shriwise, September 2019

// This is a function to test DagMC-style mesh for watertightness. For
// now this will be a stand-alone code that uses MOAB. For volumes to
// be watertight, the facet edges of each surface must be matched
// one-to-one. By default this checks for topological watertightness.
// To instead use a geometrical tolerance between two vertices that are
// considered the same, pass in a tolerance.
//
// input:  input_file h5m filename, 
//         output_file h5m filename (optional), 
//         tolerance(optional), 
//         verbose(optional)
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
#include "moab/ProgOptions.hpp"

#include "CheckWatertight.hpp"

int main(int argc, char* argv[]) {

  ProgOptions po("check_watertight: a tool for preprocessing DAGMC files to test DagMC-style mesh for watertightness");

  clock_t start_time = clock();

  bool verbose = false;
  bool check_topology;

  std::string input_file;
  std::string output_file;
  double tolerance = -1.0;

  po.addOpt<void>("verbose,v", "Verbose output", &verbose);
  
  po.addRequiredArg<std::string>("input_file", "Path to h5m DAGMC file to proccess", &input_file);
  po.addOpt<std::string>("output_file,o", "Specify the output filename (default is to overwrite in input_file)", &output_file);
  po.addOpt<double>("tolerance,t", "Specify a coincidence tolerance for triangle vertices. If no tolerance is specified, a more robust, topological check of the DAGMC mesh will occur by default.", &tolerance);

  po.parseCommandLine(argc, argv);

  if (output_file == "")
    output_file = input_file;

  static moab::Core instance;
  moab::Interface* mbi = &instance;

  // load file
  moab::ErrorCode result;
  moab::EntityHandle input_set;
  result = mbi->create_meshset(moab::MESHSET_SET, input_set);   //create handle to meshset
  if (moab::MB_SUCCESS != result) {
    return result;
  }

  result = mbi->load_file(input_file.c_str(), &input_set);   //load the file into the meshset
  if (moab::MB_SUCCESS != result) {
    // failed to load the file
    std::cout << "could not load file" << std::endl;
    return result;
  }


  if (tolerance == -1.0){
    std::cout << "geometry check" << std::endl;
    check_topology = false;
  }else if (tolerance > 0){
    std::cout << "topology check" << std::endl;
    check_topology = true;
  }else{
    MB_CHK_SET_ERR(moab::MB_FAILURE, "A proximity tolerance of " << tolerance << " was provided. Please provide a tolerance greater than or equal to zero.");
  }

  // replaced much of this code with a more modular version in check_watertight_func for testing purposes
  std::set<int> leaky_surfs, leaky_vols;
  bool sealed, test;
  test = false;
  // is the order of the optional variables going to be a problem?
  // (i.e. we 'skipped' the variable test)
  CheckWatertight cw = CheckWatertight(mbi);
  result = cw.check_mesh_for_watertightness(input_set, tolerance, sealed, test, verbose, check_topology);
  MB_CHK_SET_ERR(result, "could not check model for watertightness");

  clock_t end_time = clock();
  std::cout << (double)(end_time - start_time) / CLOCKS_PER_SEC << " seconds" << std::endl;

}


