// ********************************************************************
// Brandon Smith
// August, 2009

/* _curve_to_be_tested_for_watertightness_
      vert1 X X vert1
            | |
      vert2 X |
  surf1     | |    surf2
            | |
      vert3 X X vert2
            | |
      vert4 X X vert3                   */

// input:  h5m filename, tolerance
// output: watertight h5m

// make CXXFLAGS=-g for debug
// make CXXFLAGS=-pg for profiling

// modified by Andrew Davis 2012
// Updated deprecated MOAB calls

#include <iostream>
#include <sstream>
#include <iomanip> // for setprecision
#include <limits> // for min/max values
#include <assert.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"
#include "moab/GeomTopoTool.hpp"

#include "mw_func.hpp"
#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"
#include "cleanup.hpp"




MBInterface *MOAB();

int main(int argc, char **argv) 
  {

// ******************************************************************
    // Load the h5m file and create tags.
    // ******************************************************************

    clock_t start_time = clock();


    // check input args
    if( 2 > argc || 3 < argc ) 
      {
	std::cout << "To zip a faceted h5m file:" << std::endl;
	std::cout << "$ ./make_watertight <input_file.h5m>" << std::endl;
	std::cout << "To facet and zip an ACIS file using the default facet tolerance:" << std::endl;
	std::cout << "$ ./make_watertight <input_file.sat>" << std::endl;
	std::cout << "To facet and zip an ACIS file using a specified facet tolerance:" << std::endl;
	std::cout << "$ ./make_watertight <input_file.sat> <facet_tolerance>" << std::endl;
	return MB_FAILURE;
      }

    // The root name does not have an extension
    std::string input_name = argv[1];
    std::string root_name = argv[1];
    int len = root_name.length();
    root_name.erase(len - 4);
    bool is_acis;

    // load the input file
    MBErrorCode result, rval;
    MBEntityHandle input_set;

    rval = MBI()->create_meshset( MESHSET_SET, input_set );

    if(gen::error(MB_SUCCESS!=rval,"failed to create_meshset"))
      {
	return rval;
      }

    std::cout << "Loading input file..." << std::endl;

    // If reading an h5m file, the facet tolerance has already been determined.
    // Read the facet_tol from the file_set. There should only be one input
    // argument.

    if(std::string::npos!=input_name.find("h5m") && (2==argc)) 
      {
	rval = MBI()->load_file( input_name.c_str(), &input_set );
	if(gen::error(MB_SUCCESS!=rval,"failed to load_file 0")) 
	  {
	    return rval;      
	  }
	
	is_acis = false;

    // If reading a sat file, the facet toleance will default to 1e-3 if it is
    // not specified. If the user does not specify a facet_tol, default to 1e-3.
    // This is the same as what ReadCGM uses.
      } 

    /*
     // recreate to only perform these operations on h5m meshes  
    else if(std::string::npos!=input_name.find("sat") && 
	      ((2==argc) || (3==argc)) ) 
      {
	double facet_tol;
	if(3 == argc) 
	  {
	    facet_tol = atof(argv[2]);
	  }
	else 
	  {
	    facet_tol = 1e-3;
	  }

	std::string options;
	options += "FACET_DISTANCE_TOLERANCE=";
	std::stringstream facet_tol_ss;
	facet_tol_ss << facet_tol; 
	options += facet_tol_ss.str();
	if(debug) std::cout << "  options=" << options << std::endl;
	rval = MBI()->load_file( input_name.c_str(), &input_set, options.c_str() );
	if(gen::error(MB_SUCCESS!=rval,"failed to load_file 1")) return rval;      

      // write an HDF5 file of facets with known tolerance   
	std::string facet_tol_filename = root_name + "_" + facet_tol_ss.str() + ".h5m";
	rval = MBI()->write_mesh( facet_tol_filename.c_str() );
	if(gen::error(MB_SUCCESS!=rval,"failed to write_mesh 0")) return rval;      
	is_acis = true;
      } 
    else 
      {
	std::cout << "incorrect input arguments" << std::endl;
	return MB_FAILURE;
      }
     //not required if  only doing this with h5m files
     */
    //loading completed at this point
    clock_t load_time = clock();  

    //seal the input mesh set
    double facet_tol;
    result= mw_func::make_mesh_watertight(input_set, facet_tol);
    if(gen::error(MB_SUCCESS!=result, "could not make model watertight")) return result;
    
    //write file  
    std::string output_filename;
    if(is_acis) {  
      std::stringstream facet_tol_ss;
      facet_tol_ss << facet_tol; 
      output_filename = root_name + "_" + facet_tol_ss.str() + "_zip.h5m";
    } else {
      output_filename = root_name + "_zip.h5m";
    }
    // PROBLEM: If I write the input meshset the writer returns MB_FAILURE.
    // This happens only if I delete vertices when merging.
    // result = MBI()->write_mesh( filename_new.c_str(), &input_meshset, 1);
    clock_t zip_time = clock();
    result = MBI()->write_mesh( output_filename.c_str() );
    if (MB_SUCCESS != result) std::cout << "result= " << result << std::endl;
    assert(MB_SUCCESS == result);

    clock_t write_time = clock();
    std::cout << "Timing(seconds): loading=" 
              << (double) (load_time -start_time)/CLOCKS_PER_SEC << ", sealing=" 
              << (double) (zip_time  -load_time )/CLOCKS_PER_SEC << ", writing=" 
              << (double) (write_time-zip_time  )/CLOCKS_PER_SEC << std::endl;

  return 0;  
  }

MBInterface *MBI() 
{
    static MBCore instance;
    return &instance;
}
