// ********************************************************************
// Brandon Smith
// August, 2009

// input:  h5m filename, tolerance
// output: watertight h5m

// make CXXFLAGS=-g for debug
// make CXXFLAGS=-pg for profiling

#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

// moab includes
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include "moab/GeomTopoTool.hpp"

#include "MakeWatertight.hpp"

moab::ErrorCode write_sealed_file( moab::Interface* mbi, std::string root_filename, double facet_tol, bool is_acis);

int main(int argc, char **argv)
{

// ******************************************************************
// Load the h5m file
// ******************************************************************

  clock_t start_time = clock();

  static moab::Core instance;
  moab::Interface* mbi = &instance;

  // check input args
  if( 2 > argc || 3 < argc ) {
    std::cout << "To zip a faceted h5m file:" << std::endl;
    std::cout << "$ ./make_watertight <input_file.h5m>" << std::endl;
    std::cout << "To facet and zip an ACIS file using the default facet tolerance:" << std::endl;
    std::cout << "$ ./make_watertight <input_file.sat>" << std::endl;
    std::cout << "To facet and zip an ACIS file using a specified facet tolerance:" << std::endl;
    std::cout << "$ ./make_watertight <input_file.sat> <facet_tolerance>" << std::endl;
    return moab::MB_FAILURE;
  }

  // The root name does not have an extension
  std::string input_name = argv[1];
  std::string root_name = argv[1];
  int len = root_name.length();
  root_name.erase(len - 4);
  bool is_acis;

  // load the input file
  moab::ErrorCode result, rval;
  moab::EntityHandle input_set;

  rval = mbi->create_meshset( moab::MESHSET_SET, input_set );

  MB_CHK_SET_ERR(rval,"failed to create_meshset");

  std::cout << "Loading input file..." << std::endl;

  // If reading an h5m file, the facet tolerance has already been determined.
  // Read the facet_tol from the file_set. There should only be one input
  // argument.

  if(std::string::npos!=input_name.find("h5m") && (2==argc)) {
    rval = mbi->load_file( input_name.c_str(), &input_set );
    MB_CHK_SET_ERR(rval,"failed to load_file 0");
    is_acis = false;
  }
  //loading completed at this point
  clock_t load_time = clock();
  //seal the input mesh set
  double facet_tol;
  MakeWatertight mw(mbi);
  result= mw.make_mesh_watertight(input_set, facet_tol);
  MB_CHK_SET_ERR(result, "could not make model watertight");


  //write file
  clock_t zip_time = clock();
  std::cout << "Writing zipped file..." << std::endl;
  write_sealed_file( mbi, root_name, facet_tol, is_acis);
  MB_CHK_SET_ERR(result, "could not write the sealed mesh to a new file");


  clock_t write_time = clock();
  std::cout << "Timing(seconds): loading="
            << (double) (load_time -start_time)/CLOCKS_PER_SEC << ", sealing="
            << (double) (zip_time  -load_time )/CLOCKS_PER_SEC << ", writing="
            << (double) (write_time-zip_time  )/CLOCKS_PER_SEC << std::endl;

  return 0;
}

moab::ErrorCode write_sealed_file( moab::Interface* mbi, std::string root_filename, double facet_tol, bool is_acis)
{

  moab::ErrorCode result;
  std::string output_filename;
  if(is_acis) {
    std::stringstream facet_tol_ss;
    facet_tol_ss << facet_tol;
    output_filename = root_filename + "_" + facet_tol_ss.str() + "_zip.h5m";
  } else {
    output_filename = root_filename + "_zip.h5m";
  }
  result = mbi->write_mesh( output_filename.c_str() );
  if (moab::MB_SUCCESS != result) std::cout << "result= " << result << std::endl;
  assert(moab::MB_SUCCESS == result);

  return result;
}
