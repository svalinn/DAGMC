//
// Patrick Shriwise
// September 2013
// This program is designed to run a set of tests on the make_watertight algorithm.
// This will be a stand-alone piece of code that uses MOAB to open, modify
// (break), and re-seal geometries/ 
// input: cyl.h5m file (found in ../make_watertight/test/)
// output: pass/fail for each of the tests


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <algorithm>
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/Range.hpp"
#include "moab/Skinner.hpp"

#include "test_funcs.hpp"

#include "mw_func.hpp"
#include "cw_func.hpp"
#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"

int main(int argc, char **argv)
{

//open moab instance
moab::Interface *MBI();

//for unit testing purposes, we don't care about the output. Just PASS or FAIL. 
bool verbose=false;

// ******************************************************************
  // Load the h5m file and create tags.
  // ******************************************************************

  clock_t start_time;
  start_time = clock();
  //save the mesh to a new filename
  std::string input_name="cyl.h5m";
  std::string root_name="cyl.h5m";
  int len = root_name.length();
  root_name.erase(len-4);

  // load file and get tolerance from input argument
  moab::ErrorCode result;
  std::string filename = "cyl.h5m"; //set filename
  moab::EntityHandle input_set;
  result = MBI()->create_meshset( moab::MESHSET_SET, input_set ); //create handle to meshset
  if(moab::MB_SUCCESS != result) 
    {
      return result;
    }

result = MBI()->load_file( filename.c_str(), &input_set ); //load the file into the meshset
  if(moab::MB_SUCCESS != result) 
    {
      // failed to load the file
      std::cout << "could not load file" << std::endl;
      return result;
    }

  //// get faceting tolerance ////
  double facet_tolerance;

  moab::Tag faceting_tol_tag;
  //get faceting tolerance handle from file
  result = MBI()->tag_get_handle( "FACETING_TOL", 1, moab::MB_TYPE_DOUBLE,
        faceting_tol_tag , moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
  if(gen::error(moab::MB_SUCCESS!=result, "could not get the faceting tag handle")) return result;
  
  //get the faceting tolerance of any entity
  moab::Range file_set;
  result = MBI()->get_entities_by_type_and_tag( 0, moab::MBENTITYSET, 
                        &faceting_tol_tag, NULL, 1, file_set );

  //get facetint tolerance value
  result = MBI()->tag_get_data( faceting_tol_tag, &file_set.front(), 1, &facet_tolerance );
  if(gen::error(moab::MB_SUCCESS!=result, "could not get the faceting tolerance")) return result; 

  // create tags on geometry
  moab::Tag geom_tag, id_tag;
  result = MBI()->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, 
                            moab::MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  if(gen::error(moab::MB_SUCCESS != result, "could not get GEOM_DIMENSION_TAG_NAME handle")) return result; 

  result = MBI()->tag_get_handle( GLOBAL_ID_TAG_NAME, 1, 
                            moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  if(gen::error(moab::MB_SUCCESS != result, "could not get GLOBAL_ID_TAG_NAME handle")) return result;
  
  // get surface and volume sets
  moab::Range surf_sets, vol_sets; // moab::Range of set of surfaces and volumes
  // surface sets
  int dim = 2;
  void* input_dim[] = {&dim};
  result = MBI()->get_entities_by_type_and_tag( input_set, moab::MBENTITYSET, &geom_tag, 
                                                input_dim, 1, surf_sets);
  if(moab::MB_SUCCESS != result) 
    {
      return result;
    }

  // volume sets
  dim = 3;
  result = MBI()->get_entities_by_type_and_tag( input_set, moab::MBENTITYSET, &geom_tag, 
                                                input_dim, 1, vol_sets);
  if(moab::MB_SUCCESS != result)
    {
      return result;
    }
    
  //vertex sets
  dim= 0;
  moab::Range verts;
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;
  
  if(gen::error(moab::MB_SUCCESS!=result, "could not get vertex coordinates")) return result; 
  
if(verbose)
{
  std::cout<< "number of verticies= " << verts.size() << std::endl;  
  std::cout<< "number of surfaces= " << surf_sets.size() << std::endl;
  std::cout<< "number of volumes= "  << vol_sets.size() << std::endl;
}

//initialize booleans to pass to make_mesh_watertight

  bool check_topology, test, sealed;
       check_topology=false;
       test=true;
  
// initialize boolean for each set of tests
  bool test_set_result=true;

///////////Single Verticie Movement Tests////////////////////
  std::cout << "SINGLE VERTEXT MOVEMENT TESTS" << std::endl;

///////////////BEGIN 1st TEST////////////////////////

 std::cout << "Test 1: vertex bump in x direction:";

  //perform single vertex bump and test
  result=single_vert_bump(verts, 0.9*facet_tolerance , 0 , 0 ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // seal the model using make_watertight 
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;
  
  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
  if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 2ND TEST////////////////////////

  std::cout << "Test 2: vertex bump in y direction:";

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  // "Break" the geometry
  result=single_vert_bump(verts, 0 , 0.9*facet_tolerance , 0); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

// Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }


///////////////BEGIN 3RD TEST////////////////////////

  std::cout << "Test 3: vertex bump in z direction:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  // "Break" the geometry
  result = single_vert_bump(verts, 0 , 0 , 0.9*facet_tolerance ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // Seal the mesh
  result = mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

// Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 4TH TEST////////////////////////

  std::cout << "Test 4: vertex bump in rand direction:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=rand_vert_bump(verts, facet_tolerance, root_name); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create single vertex bump test")) return result;
  
  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

// Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

  if(sealed)
  {
   std::cout << "PASS" << std::endl << std::endl ;
  }
  else
  {
   std::cout << "FAIL" << std::endl << std::endl ;
   test_set_result=false;
  }

//////////END SINGLE VERTEX MOVEMENT TESTS///////////

  if(test_set_result)
  {
   std::cout << "SINGLE VERTEX MOVEMENT TESTS PASSED" << std::endl << std::endl;
  }
  else
  {
   std::cout << "SINGLE VERTEX MOVEMENT TESTS FAILED" << std::endl << std::endl;
   exit(0);
  }

//////////LOCKED VERTEX PAIR MOVEMENT TESTS/////////////////

std::cout << "LOCKED VERTEX PAIR MOVEMENT TESTS" << std::endl;

///////////////BEGIN 1ST TEST////////////////////////

  std::cout << "Test 1: locked pair of verticies move in x:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=locked_pair_bump(verts, 0.9*facet_tolerance , 0 , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;


   // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 2ND TEST////////////////////////

  std::cout << "Test 2: locked pair of verticies move in y:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=locked_pair_bump(verts, 0 , 0.9*facet_tolerance , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;


   // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }


///////////////BEGIN 3RD TEST////////////////////////

  std::cout << "Test 3: locked pair of verticies move in z:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=locked_pair_bump(verts, 0 , 0 , 0.9*facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;


   // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 4TH TEST////////////////////////

  std::cout << "Test 4: adj. pair of verticies random move:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=locked_pair_bump_rand(verts, facet_tolerance , root_name); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  

   // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;


  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  

  if(sealed)
  {
   std::cout << "PASS" << std::endl << std::endl ;
  }
  else
  {
   std::cout << "FAIL" << std::endl << std::endl ;
   test_set_result=false;
  }


//////////END LOCKED VERTEX PAIR MOVEMENT TESTS///////////

  if(test_set_result)
  {
   std::cout << "LOCKED VERTEX PAIR MOVEMENT TESTS PASSED" << std::endl << std::endl;
  }
  else
  {
   std::cout << "LOCKED VERTEX PAIR MOVEMENT TESTS FAILED" << std::endl << std::endl;
   exit(0);
  }


////////////RAND PAIR MOVEMENT TESTS//////////////////

std::cout << "RAND PAIR MOVEMENT TESTS" << std::endl;

///////////////BEGIN 1ST TEST////////////////////////

  std::cout << "Test 1: random locked pair of verticies move in x:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;
  
  // "Break" the geometry
  result=rand_locked_pair_bump(verts, 0.9*facet_tolerance , 0 , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 2ND TEST////////////////////////

  std::cout << "Test 2: random locked pair of verticies move in y:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  
  // "Break" the geometry
  result=rand_locked_pair_bump(verts, 0 , 0.9*facet_tolerance , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 3RD TEST////////////////////////

  std::cout << "Test 3: random locked pair of verticies move in z:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;
  
  // "Break" the geometry
  result=rand_locked_pair_bump(verts, 0 , 0 , 0.9*facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }


///////////////BEGIN 4TH TEST////////////////////////

  std::cout << "Test 4: random pair of verticies move in random dir:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  // "Break" the geometry
  result=rand_locked_pair_bump_rand(verts, facet_tolerance , root_name); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl << std::endl ;
  }
  else
  {
   std::cout << "FAIL" << std::endl << std::endl ;
   test_set_result=false;
  }


//////////END RAND LOCKED VERTEX PAIR MOVEMENT TESTS///////////

  if(test_set_result)
  {
   std::cout << "RANDOM LOCKED VERTEX PAIR MOVEMENT TESTS PASSED" << std::endl << std::endl;
  }
  else
  {
   std::cout << "RANDOM LOCKED VERTEX PAIR MOVEMENT TESTS FAILED" << std::endl << std::endl;
   exit(0);
  }


/////////////ADJACENT PLUS ONE TESTS//////////////////

std::cout << "ADJACENT PLUS ONE TESTS" << std::endl;

///////////////BEGIN 1ST TEST////////////////////////

  std::cout << "Test 1: locked pair of verticies (adj+1) move in x:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result=adjplone_locked_pair_bump(verts, 0.9*facet_tolerance , 0 , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 2ND TEST////////////////////////

  std::cout << "Test 2: locked pair of verticies (adj+1) move in y:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result=adjplone_locked_pair_bump(verts, 0 , 0.9*facet_tolerance , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 3RD TEST////////////////////////

  std::cout << "Test 3: locked pair of verticies (adj+1) move in z:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result=adjplone_locked_pair_bump(verts, 0 , 0 , 0.9*facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }


///////////////BEGIN 4TH TEST////////////////////////

  std::cout << "Test 4: pair of verticies (adj+1) move in random dir:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result=adjplone_locked_pair_bump_rand(verts, facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl << std::endl ;
  }
  else
   
  {
   std::cout << "FAIL" << std::endl << std::endl ;
   test_set_result=false;
  }

//////////END ADJACENT + 1 VERTEX MOVEMENT TESTS///////////

  if(test_set_result)
  {
   std::cout << "ADJACENT PLUS ONE VERTEX PAIR MOVEMENT TESTS PASSED" << std::endl << std::endl;
  }
  else
  {
   std::cout << "ADJACENT PLUS ONE VERTEX PAIR MOVEMENT TESTS FAILED" << std::endl << std::endl;
   exit(0);
  }


//////////NON-ADJACENT LOCKED PAIR TESTS//////////////

std::cout << "NON-ADJACENT LOCKED PAIR TESTS" << std::endl;

///////////////BEGIN 1ST TEST////////////////////////

  std::cout << "Test 1: non-adjacent locked pair of verticies move in x:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result= nonadj_locked_pair_bump(verts,  0.9*facet_tolerance , 0 , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 2ND TEST////////////////////////

  std::cout << "Test 2: non-adjacent locked pair of verticies move in y:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result= nonadj_locked_pair_bump(verts, 0 ,  0.9*facet_tolerance , 0 , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }


///////////////BEGIN 3RD TEST////////////////////////

  std::cout << "Test 3: non-adjacent locked pair of verticies move in z:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result= nonadj_locked_pair_bump(verts, 0 , 0 ,  0.9*facet_tolerance  , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl;
  }
  else
  {
   std::cout << "FAIL" << std::endl;
   test_set_result=false;
  }

///////////////BEGIN 4TH TEST////////////////////////

  std::cout << "Test 4: non-adjacent pair of verticies move in random dir:" ;

  // Clear the mesh and reload original geometry for the next test
  result=reload_mesh( filename.c_str(), input_set);
  if (gen::error(moab::MB_SUCCESS!=result, "could not reload the mesh" )) return result; 

  // retrieve the verticies again so the model can be broken
  result = MBI()->get_entities_by_dimension(input_set, dim, verts, false);
  if(gen::error(moab::MB_SUCCESS!=result, " could not get vertices from the mesh")) return result;

  //"Break" the geometry
  result= nonadj_locked_pair_bump_rand(verts, facet_tolerance , root_name ); 
  if(gen::error(moab::MB_SUCCESS!=result, "could not create locked pair vertex bump test")) return result;

  // Seal the mesh
  result=mw_func::make_mesh_watertight (input_set, facet_tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not make the mesh watertight")) return result;

  // Lastly Check to see if make_watertight fixed the model
  result=cw_func::check_mesh_for_watertightness( input_set, facet_tolerance, sealed, test);
  if(gen::error(moab::MB_SUCCESS!=result, "could not check model for watertightness")) return result;
  
    if(sealed)
  {
   std::cout << "PASS" << std::endl << std::endl ;
  }
  else
  {
   std::cout << "FAIL" << std::endl << std::endl ;
   test_set_result=false;
  }


//////////END NON-ADJACNET LOCKED VERTEX PAIR MOVEMENT TESTS///////////

  if(test_set_result)
  {
   std::cout << "NON-ADJACENT LOCKED VERTEX PAIR MOVEMENT TESTS PASSED" << std::endl << std::endl;
  }
  else
  {
   std::cout << "NON-ADJACENT LOCKED VERTEX PAIR MOVEMENT TESTS FAILED" << std::endl << std::endl;
   exit(0);
  }


}


