///
/// Patrick Shriwise
/// September 2013
/// This program is designed to run a set of tests on the make_watertight algorithm.
/// This will be a stand-alone piece of code that uses MOAB to open, modify
/// (break), and re-seal geometries/ 
/// input: cyl.h5m file (found in ../make_watertight/test/)
/// output: pass/fail for each of the tests


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

#include "mw_func.hpp"
#include "cw_func.hpp"
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

MBErrorCode reload_mesh( MBEntityHandle meshset, const char* filename, bool debug = false) {

MBErrorCode result;

  // clear meshset
  result= MBI() -> clear_meshset(&meshset,1);
  if(gen::error(MB_SUCCESS!=result, "could not clear the mesh set")) return result;

  // ensure that the meshset no longer contains entities
  int num_meshsets;
  result= MBI() -> num_contained_meshsets ( meshset, &num_meshsets);
  if(gen::error(MB_SUCCESS!=result, "could not get the number of meshsets")) return result; 

  if (debug) std::cout << "number of mesh sets after clear = " << num_meshsets << std::endl;
  if (num_meshsets!=0) return MB_FAILURE;
  
  // re-initialize meshset
  result= MBI() -> create_meshset(MESHSET_SET, meshset);
  if(gen::error(MB_SUCCESS!=result, "could not create the new meshset")) return result; 
 
  //reload the file
  result = MBI() -> load_file(filename, &meshset);
  if(gen::error(MB_SUCCESS!=result, "could not re-load the file into the mesh set")) return result; 
   
  //check that something was loaded into the meshset
  result= MBI() -> num_contained_meshsets ( meshset, &num_meshsets);
  if(gen::error(MB_SUCCESS!=result, "could not get the number of meshsets")) return result; 

  if(debug) std::cout << "number of mesh sets after reload = " << num_meshsets << std::endl;

  if(num_meshsets==0) return MB_FAILURE;

 
  

return result;
}



MBErrorCode single_vert_bump( MBRange verts, double bump_dist_x, double bump_dist_y, double bump_dist_z,  std::string root_name, bool verbose = true ) {
 
  MBErrorCode result; 
  
  //get vertex coordinates
  double coords[3];
  MBEntityHandle vertex=verts.back();
  result= MBI()-> get_coords( &vertex , 1 , coords );
 
 if(verbose)
 {
  //original coordinates
  std::cout << std::endl << "Original Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;

  coords[0]+=bump_dist_x;
  coords[1]+=bump_dist_y;
  coords[2]+=bump_dist_z;
  
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
 }
  //write new coordinates to the mesh
  // might not be necesarry any longer as we move to doing tests on a moab-instance basis
  result=MBI()-> set_coords( &vertex, 1, coords);
  if (gen::error(MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;
  // alter output filename
  std::string output_filename = root_name + "_mod.h5m";
  //write file
  result=MBI()->write_mesh(output_filename.c_str());
  if(gen::error(MB_SUCCESS!=result,"could not write the mesh to output_filename")) return result; 

  return MB_SUCCESS;
}

int main(int argc, char **argv)
{

//open moab instance
MBInterface *MBI();

//for unit testing purposes, we don't care about the output. Just PASS or FAIL. 
bool verbose=false;

// ******************************************************************
  // Load the h5m file and create tags.
  // ******************************************************************

  clock_t start_time;
  start_time = clock();
  //save the mesh to a new filename
  std::string input_name=argv[1];
  std::string root_name=argv[1];
  int len = root_name.length();
  root_name.erase(len-4);
  // check input args
  if( argc < 2 || argc > 5) 
    {
    std::cout << "To check using topology of facet points:              " << std::endl;
    std::cout << "./test_cyl <filename> <verbose(true or false)>" << std::endl;
    std::cout << "To check using geometry tolerance of facet points:    " << std::endl;
    std::cout << "./test_cyl <filename> <verbose(true or false)> <tolerance>" << std::endl;
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

  //standard tolerance for tests done on cylinders
  double tol=1e-04;

  // create tags on geometry
  MBTag geom_tag, id_tag;
  result = MBI()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, 
                            MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  if(gen::error(MB_SUCCESS != result, "could not get GEOM_DIMENSION_TAG_NAME handle")) return result; 

  result = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1, 
                            MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  if(gen::error(MB_SUCCESS != result, "could not get GLOBAL_ID_TAG_NAME handle")) return result;

  
  // get surface and volume sets
  MBRange surf_sets, vol_sets; // MBRange of set of surfaces and volumes
  // surface sets
  int dim = 2;
  void* input_dim[] = {&dim};
  result = MBI()->get_entities_by_type_and_tag( input_set, MBENTITYSET, &geom_tag, 
                                                input_dim, 1, surf_sets);
  if(MB_SUCCESS != result) 
    {
      return result;
    }

  // volume sets
  dim = 3;
  result = MBI()->get_entities_by_type_and_tag( input_set, MBENTITYSET, &geom_tag, 
                                                input_dim, 1, vol_sets);
  if(MB_SUCCESS != result)
    {
      return result;
    }
    
  //vertex sets
  dim= 0;
  MBRange verts;
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;
  
  if(gen::error(MB_SUCCESS!=result, "could not get vertex coordinates")) return result; 

if(verbose)
{
  std::cout<< "number of verticies= " << verts.size() << std::endl;  
  std::cout<< "number of surfaces= " << surf_sets.size() << std::endl;
  std::cout<< "number of volumes= "  << vol_sets.size() << std::endl;
}

  //perform single vertex bump and test
  result=single_vert_bump(verts, 0 , 0 , 10*tol, root_name, verbose ); 
  if(gen::error(MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
 
  
  bool check_topology, test, sealed;
  check_topology=false;
  test=true;
  
  /*
  //check that geometry is unsealed
  result=cw_func::check_mesh_for_watertightness( input_set, tol, sealed, test);
  if(gen::error(MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  */

  if(sealed)
  {
   std::cout << "Warning: geometry was not broken by single bump vert" << std::endl;
  }

  //seal the mesh
  double facet_tol;
  result=mw_func::make_mesh_watertight (input_set, facet_tol, verbose);
  if(gen::error(MB_SUCCESS!=result, "could not make the mesh watertight")) return result;
  
  




  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, tol, sealed, test);
  if(gen::error(MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
  if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   exit(0);
  }
  
  result=reload_mesh(input_set, filename.c_str(), true);
  if(gen::error(MB_SUCCESS!=result, "could not reload the mesh")) return result; 

 

}

MBInterface* MBI() {
 static MBCore instance;
 return &instance;
}

