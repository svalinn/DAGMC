//
// Patrick Shriwise
// September 2013
// This program is designed to run a set of tests on the make_watertight algorithm.
// This will be a stand-alone piece of code that uses MOAB to open, modify
// (break), and re-seal geometries
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

#include "gtest/gtest.h"
#include "test_funcs.hpp"

/// moves a vertex along the rim of the cylinder in the theta direction a distance equal to the faceting_tolerance
moab::ErrorCode move_vert_theta(moab::EntityHandle vertex, double tolerance, bool verbose = false) {

  moab::ErrorCode result; 
  
  //get vertex coordinates
  double coords[3];
  result = MBI()->get_coords(&vertex, 1, coords);
 
  // determine radius
  double radius = sqrt(pow(coords[0],2)+pow(coords[1],2));

  // get the current theta value
  // need both because of the oddness/evenness of the inverse functions
  double theta_x = acos(coords[0]/radius);
  double theta_y = asin(coords[1]/radius);
  
  // set the vertex bump distance
  double dtheta = tolerance/(radius);

  if(verbose) {
    //original coordinates
    std::cout << std::endl << "Original Coordinates" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }
  // create new x and y values
  coords[0] = radius*cos(theta_x+dtheta); 
  coords[1] = radius*sin(theta_y+dtheta);

  if(verbose) {
    //altered coordinates
    std::cout << std::endl << "Modified Coordinates" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }
  //set new vertex coordinates  
  result = MBI()->set_coords(&vertex, 1, coords);
  if (gen::error(moab::MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;

  return moab::MB_SUCCESS;
}

/// moves the vertex in R some distance less than tol
moab::ErrorCode move_vert_R(moab::EntityHandle vertex, double tol, bool verbose = false) {

  moab::ErrorCode result; 
  
  //get vertex coordinates
  double coords[3];
  result = MBI()->get_coords(&vertex, 1, coords);
 
  if(verbose){
  //original coordinates
  std::cout << "Vertex ID: " << gen::geom_id_by_handle(vertex) << std::endl;
  std::cout << std::endl << "Original Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
  }
  double radius = sqrt(pow(coords[0],2)+pow(coords[1],2));
  //get unit vector in x-y plane
  coords[0]/=radius;
  coords[1]/=radius;
  
  //alter radius to new value of radius+tol
  radius-=tol;
  coords[0]*=radius;
  coords[1]*=radius;
 
  if(verbose){
  //altered coordinates
  std::cout << std::endl << "Modified Coordinates" << std::endl;
  std::cout << "x = " << coords[0] << std::endl;
  std::cout << "y = " << coords[1] << std::endl;
  std::cout << "z = " << coords[2] << std::endl;
  }
  
  //set new vertex coordinates  
  result = MBI()->set_coords(&vertex, 1, coords);
  if (gen::error(moab::MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;

return moab::MB_SUCCESS;
}

/// bumps the last vertex in the cylinder model in the R direction
moab::ErrorCode single_vert_bump_R(moab::Range verts, double facet_tol, bool verbose = false) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex=verts.back();
  //move the desired vertex by the allotted distance
  result = move_vert_R(vertex, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move single vert")) return result;

  return moab::MB_SUCCESS;
}


/// selects a random pair of verticies and moves them along theta a distance less than the faceting tolerance
moab::ErrorCode rand_locked_pair_bump_theta(moab::Range verts, double facet_tol, bool verbose = false) {
 
  moab::ErrorCode result; 

  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
 
  int index = 1+static_cast<int>(num*(number_of_verts-1));

  moab::EntityHandle vertex1=(verts.back()-index);
  moab::EntityHandle vertex2=(verts.back()-index-1);
  
  //move the desired verticies by the allotted distance(s)
  result = move_vert_theta(vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_theta(vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

// FOR CYLINDER TESTING ONLY 
/// moves the last vertex in the model along the curve of the cylinder some distance bump distance theta
moab::ErrorCode theta_vert_bump(moab::Range verts, double bump_dist_theta, bool verbose = false) {
 
  moab::ErrorCode result; 
  
  //get vertex coordinates
  double coords[3];
  moab::EntityHandle vertex = verts.back();
  result = MBI()->get_coords(&vertex, 1, coords);
 
  // determine radius
  double radius = sqrt(pow(coords[0],2)+pow(coords[1],2));

  // get the current theta value
  double theta = asin(coords[1]/radius);
  
  // set the vertex bump distance
  double dtheta = 0.5*bump_dist_theta/(radius);
 
  if(verbose) {
    //original coordinates
    std::cout << std::endl << "Original Coordinates HEre I am" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }
  // create new x and y values
  coords[0] = radius*cos(theta+dtheta); 
  coords[1] = radius*sin(theta+dtheta);
  if(verbose) {
    //altered coordinates
    std::cout << std::endl << "Modified Coordinates" << std::endl;
    std::cout << "x = " << coords[0] << std::endl;
    std::cout << "y = " << coords[1] << std::endl;
    std::cout << "z = " << coords[2] << std::endl;
  }
  //write new coordinates to the mesh
  // might not be necesarry any longer as we move to doing tests on a moab-instance basis
  result = MBI()->set_coords(&vertex, 1, coords);
  if (gen::error(moab::MB_SUCCESS!=result, "could not set the vertex coordinates")) return result;
  // alter output filename
 
  return moab::MB_SUCCESS;
}

/// moves two adjacent vertices along theta a distance equal to the faceting tolerance
moab::ErrorCode locked_pair_bump_theta(moab::Range verts, double tolerance, bool verbose = false) {

  moab::ErrorCode result;

  //get vertex coordinates
  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = verts.back()-1;
  
  result = move_vert_theta(vertex1, tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result,"could not move vertex1 along theta")) return result; 
  result = move_vert_theta(vertex2, tolerance, verbose);
  if(gen::error(moab::MB_SUCCESS!=result,"could not move vertex1 along theta")) return result; 

  return moab::MB_SUCCESS;
}


/// moves the third to last and the last verticies in the model in theta the same distance along theta equal to the faceting tolerance
moab::ErrorCode adjplone_locked_pair_bump_theta(moab::Range verts, double facet_tol, bool verbose = false) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-2);
 
  //move the desired verticies by the allotted distance(s)
  result = move_vert_theta(vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_theta(vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

moab::ErrorCode locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose = false) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-1);
 
  //move the desired verticies by the allotted distance(s)
  result = move_vert_R(vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_R(vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

/// selects random verticies from verts and moves them in R a distance equal to the faceting tolerance
moab::ErrorCode rand_locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose = false) {
 
  moab::ErrorCode result; 
  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
 
  int index = 1+static_cast<int>(num*(number_of_verts-1));

  moab::EntityHandle vertex1=(verts.back()-index);
  moab::EntityHandle vertex2=(verts.back()-index-1);  
  //move the desired verticies by the allotted distance(s)
  result = move_vert_R(vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_R(vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

/// selects a the last vertex and third to last vertex in the model and moves them in R a distance equal to the faceting tolerance
moab::ErrorCode adjplone_locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose = false) {
 
  moab::ErrorCode result; 
  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-2);
 
  //move the desired verticies by the allotted distance(s)
  result = move_vert_R(vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_R(vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

/// moves the last vertex in the model and a randomly selected, non-adjacent vertex and moves them both in R a distance equal to the faceting tolerance
moab::ErrorCode nonadj_locked_pair_bump_R(moab::Range verts, double facet_tol, bool verbose = false) {
 
  moab::ErrorCode result; 
  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  int index = static_cast<int>(num*((number_of_verts-2)));

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-index);
  //move the desired verticies by the allotted distance(s)
  result = move_vert_R(vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_R(vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}

/// selects a random pair of adjacent verticies and bumps them along the theta direction a distance equal to the faceting tolerance
moab::ErrorCode nonadj_locked_pair_bump_theta(moab::Range verts, double facet_tol, bool verbose = false) {
 
  moab::ErrorCode result; 
  //select random verticies from verts
  int number_of_verts = verts.size();
  double num = rand()/static_cast<double>(RAND_MAX);
  int index = static_cast<int>(num*((number_of_verts-2)));

  moab::EntityHandle vertex1 = verts.back();
  moab::EntityHandle vertex2 = (verts.back()-index);
  
  //move the desired verticies by the allotted distance(s)
  result = move_vert_theta(vertex1, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex1")) return result;
  result = move_vert_theta(vertex2, facet_tol, verbose);
  if(gen::error(moab::MB_SUCCESS!=result, "could not move vertex2")) return result;

  return moab::MB_SUCCESS;
}


//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class MakeWatertightCylinderTest : public ::testing::Test
{
protected:
  virtual void SetUp() {
    filename = "cyl.h5m";
    reload_mesh();
  };

  void reload_mesh() {
    // delete meshset
    result = MBI()->delete_mesh();
    EXPECT_EQ(result,moab::MB_SUCCESS);

    // re-initialize meshset
    result = MBI()->create_meshset(moab::MESHSET_SET, input_fileset);
    EXPECT_EQ(result,moab::MB_SUCCESS);
 
    //reload the file
    result = MBI()->load_file(filename.c_str(), &input_fileset);
    EXPECT_EQ(result,moab::MB_SUCCESS);

    //// get faceting tolerance //// 
    moab::Tag faceting_tol_tag;
    //get faceting tolerance handle from file
    result = MBI()->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE,
				   faceting_tol_tag , moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
    EXPECT_EQ(result,moab::MB_SUCCESS);
  
    //get the faceting tolerance of any entity
    moab::Range file_set;
    result = MBI()->get_entities_by_type_and_tag(0, moab::MBENTITYSET, 
						 &faceting_tol_tag, NULL, 1, file_set);

    //get facetint tolerance value
    result = MBI()->tag_get_data(faceting_tol_tag, &file_set.front(), 1, &facet_tol);
    EXPECT_EQ(result,moab::MB_SUCCESS);
  
    //check that something was loaded into the meshset
    int num_meshsets;   
    result = MBI()->num_contained_meshsets (input_fileset, &num_meshsets);
    EXPECT_EQ(result,moab::MB_SUCCESS);  
    if(num_meshsets == 0) MB_CHK_ERR_CONT(moab::MB_FAILURE);

    //retrieve the verticies again so the model can be broken
    int dim = 0;
    result = MBI()->get_entities_by_dimension(input_fileset, dim, verts, false);
    EXPECT_EQ(result,moab::MB_SUCCESS);
  };
  
  virtual void TearDown() {
    result = MBI()->delete_mesh();
    EXPECT_EQ(result,moab::MB_SUCCESS);
  };

protected:
  std::string filename;
  moab::ErrorCode result;
  moab::EntityHandle input_fileset;
  moab::Range verts;
  double facet_tol;
};

TEST_F(MakeWatertightCylinderTest, SingleVertexMoveTests)
{
  //bump vert in x direction
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.9*facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump vert in y direction
  reload_mesh();
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.0, 0.9*facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump vert in z direction
  reload_mesh();
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.0, 0.0, 0.9*facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump vert away from the center of the cylinder
  reload_mesh();
  EXPECT_NO_THROW(result = single_vert_bump_R(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump vert along theta of the cylinder
  reload_mesh();
  EXPECT_NO_THROW(result = theta_vert_bump(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump vert in a random direction
  reload_mesh();
  EXPECT_NO_THROW(result = rand_vert_bump(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}


TEST_F(MakeWatertightCylinderTest, LockedVertexPairMoveTests)
{
  //bump a pair of verts in the same amount in the x direction
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.9*facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump a pair of verts in the same amount in the y direction
  reload_mesh();
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.0, 0.9*facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump a pair of verts in the same amount in the z direction
  reload_mesh();
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.0, 0.0, 0.9*facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump a pair of verts in the same amount in the radial direction
  reload_mesh();
  EXPECT_NO_THROW(result = locked_pair_bump_R(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump a pair of verts in the same amount in the theta direction
  reload_mesh();
  EXPECT_NO_THROW(result = locked_pair_bump_theta(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump a pair of verts in the same amount in a random direction  
  reload_mesh();
  EXPECT_NO_THROW(result = locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));

}

TEST_F(MakeWatertightCylinderTest, RandomVertexPairMoveTests)
{
  //move a randomly selected locked pair of vertices in the x direction
  EXPECT_NO_THROW(result = rand_locked_pair_bump(verts, 0.9*facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //move a randomly selected locked pair of vertices in the y direction
  reload_mesh();
  EXPECT_NO_THROW(result = rand_locked_pair_bump(verts, 0.0, 0.9*facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //move a randomly selected locked pair of vertices in the z direction  
  reload_mesh();
  EXPECT_NO_THROW(result = rand_locked_pair_bump(verts, 0.0, 0.0, 0.9*facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //move a randomly selected locked pair of vertices in the radial direction
  reload_mesh();
  EXPECT_NO_THROW(result = rand_locked_pair_bump_R(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));  
  //move a randomly selected locked pair of vertices along the theta direction
  reload_mesh();
  EXPECT_NO_THROW(result = rand_locked_pair_bump_theta(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));  
  //move a randomly selected locked pair of vertices in a random direction  
  reload_mesh();
  EXPECT_NO_THROW(result = rand_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  
}

TEST_F(MakeWatertightCylinderTest, AdjacentPlusOneVertexPairTests)
{
  //bump a vertex pair with a single, stationary vertex between them in the x direction 
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump(verts, 0.9*facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump a vertex pair with a single, stationary vertex between them in the y direction 
  reload_mesh();
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump(verts, 0.0, 0.9*facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump a vertex pair with a single, stationary vertex between them in the z direction 
  reload_mesh();
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump(verts, 0.0, 0.0, 0.9*facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump a vertex pair with a single, stationary vertex between them in the radial direction 
  reload_mesh();
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump_R(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));  
  //bump a vertex pair with a single, stationary vertex between them along the theta direction
  reload_mesh();
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump_theta(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));    
  //bump a vertex pair with a single, stationary vertex between them in a random direction  
  reload_mesh();
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}

TEST_F(MakeWatertightCylinderTest, NonAdjacentVertexPairTests)
{
  //move two randomly selected, non-adjacent vertices in the x direction
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump(verts, 0.9*facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //move two randomly selected, non-adjacent vertices in the y direction
  reload_mesh();
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump(verts, 0.0, 0.9*facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //move two randomly selected, non-adjacent vertices in the z direction
  reload_mesh();
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump(verts, 0.0, 0.0, 0.9*facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //move two randomly selected, non-adjacent vertices in the radial direction
  reload_mesh();
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump_R(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //move two randomly selected, non-adjacent vertices along the theta direction
  reload_mesh();
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump_theta(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //move two randomly selected, non-adjacent vertices in a random direction  
  reload_mesh();
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
