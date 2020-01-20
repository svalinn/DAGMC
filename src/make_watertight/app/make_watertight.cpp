// ********************************************************************
// input:  input_file h5m filename,
//         output_file h5m filename (optional),
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
#include "moab/ProgOptions.hpp"

#include "MakeWatertight.hpp"
#include "CheckWatertight.hpp"

moab::ErrorCode write_sealed_file(moab::Interface* mbi, std::string output_file);

int main(int argc, char* argv[]) {

  ProgOptions po("make_watertight: a tool for preprocessing DAGMC files to seal a faceted h5m file");

  clock_t start_time = clock();

  std::string input_file;
  std::string output_file;

  po.addRequiredArg<std::string>("input_file", "Path to h5m DAGMC file to proccess", &input_file);
  po.addOpt<std::string>("output_file,o", "Specify the output filename (default watertight_dagmc.h5m)", &output_file);

  po.parseCommandLine(argc, argv);

  if (output_file == "")
    output_file = "watertight_dagmc.h5m";

  static moab::Core instance;
  moab::Interface* mbi = &instance;

  // load the input file
  moab::ErrorCode result, rval;
  moab::EntityHandle input_set;

  rval = mbi->create_meshset(moab::MESHSET_SET, input_set);

  MB_CHK_SET_ERR(rval, "failed to create_meshset");

  std::cout << "Loading input file..." << std::endl;

  // When reading an h5m file, the facet tolerance has already been determined.
  // Read the facet_tol from the file_set.

  rval = mbi->load_file(input_file.c_str(), &input_set);
  MB_CHK_SET_ERR(rval, "Failed to open file: " << input_file);

  //loading completed at this point
  clock_t load_time = clock();
  //seal the input mesh set
  double facet_tol;
  MakeWatertight mw(mbi);
  result = mw.make_mesh_watertight(input_set, facet_tol);
  MB_CHK_SET_ERR(result, "could not make model watertight");


  //write file
  clock_t seal_time = clock();
  std::cout << "Writing sealed file..." << std::endl;
  write_sealed_file(mbi, output_file);
  MB_CHK_SET_ERR(result, "could not write the sealed mesh to a new file");

  clock_t write_time = clock();

  bool sealed;
  CheckWatertight cw = CheckWatertight(mbi);
  result = cw.check_mesh_for_watertightness(input_set, -1.0, sealed);
  MB_CHK_SET_ERR(result, "could not check model for watertightness");

  std::cout << "Topological check performed" << std::endl;
  std::cout << "To perform a tolerance based check run the check_Watertight tool" << std::endl;
  std::cout << "$ check_watertight <dagmc_file> -t <tolerance> " << std::endl;


  clock_t check_time = clock();
  std::cout << "Timing(seconds): loading="
            << (double)(load_time - start_time) / CLOCKS_PER_SEC << ", sealing="
            << (double)(seal_time  - load_time) / CLOCKS_PER_SEC << ", writing="
            << (double)(write_time - seal_time) / CLOCKS_PER_SEC << ", checking="
            << (double)(check_time - write_time) / CLOCKS_PER_SEC << std::endl;

  return 0;
}

moab::ErrorCode write_sealed_file(moab::Interface* mbi, std::string output_file) {

  moab::ErrorCode result;

  result = mbi->write_mesh(output_file.c_str());
  if (moab::MB_SUCCESS != result)
    std::cout << "result= " << result << std::endl;
  assert(moab::MB_SUCCESS == result);

  return result;
}
