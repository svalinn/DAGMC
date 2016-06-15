//
// Patrick Shriwise
// June 2016
// This program is designed to run a set of tests on the make_watertight algorithm.
// This will be a stand-alone piece of code that uses MOAB to open, modify
// (break), and re-seal geometries
// input: cones.h5m
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

#include "gen.hpp"
#include "arc.hpp"
#include "zip.hpp"

#include "gtest/gtest.h"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class MakeWatertightConeTest : public ::testing::Test
{
protected:
  virtual void SetUp() {
    filename = "cones.h5m";
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

TEST_F(MakeWatertightConeTest, SingleVertexMoveTests)
{
  //bump in x direction
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.9*facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump in y direction
  reload_mesh();
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.0, 0.9*facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump in z direction
  reload_mesh();
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.0, 0.0, 0.9*facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump in random direction
  reload_mesh();
  EXPECT_NO_THROW(result = rand_vert_bump(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol)); 
}

TEST_F(MakeWatertightConeTest, LockedVertexPairMoveTests)
{
  //bump locked pair in x direction
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.9*facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump locked pair in y direction  
  reload_mesh();
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.0, 0.9*facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump locked pair in z direction
  reload_mesh();
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.0, 0.0, 0.9*facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump locked pair in random diretion
  reload_mesh();
  EXPECT_NO_THROW(result = locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));

}

TEST_F(MakeWatertightConeTest, RandomVertexPairMoveTests)
{
  //bump random locked pair in x direction
  EXPECT_NO_THROW(result = rand_locked_pair_bump(verts, 0.9*facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump random locked pair in y direction
  reload_mesh();
  EXPECT_NO_THROW(result = rand_locked_pair_bump(verts, 0.0, 0.9*facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump random locked pair in z direction  
  reload_mesh();
  EXPECT_NO_THROW(result = rand_locked_pair_bump(verts, 0.0, 0.0, 0.9*facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump random locked pair in random direction
  reload_mesh();
  EXPECT_NO_THROW(result = rand_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}

TEST_F(MakeWatertightConeTest, AdjacentPlusOneVertexPairTests)
{
  //bump pair of verices in the x direction
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump(verts, 0.9*facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump pair of verices in the y direction
  reload_mesh();
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump(verts, 0.0, 0.9*facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump pair of verices in the z direction
  reload_mesh();
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump(verts, 0.0, 0.0, 0.9*facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump pair of verices in a random direction
  reload_mesh();
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}

TEST_F(MakeWatertightConeTest, NonAdjacentVertexPairTests)
{
  //bump a non-adjacent pair of vertices in the x direction
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump(verts, 0.9*facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump a non-adjacent pair of vertices in the y direction
  reload_mesh();
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump(verts, 0.0, 0.9*facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump a non-adjacent pair of vertices in the z direction
  reload_mesh();
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump(verts, 0.0, 0.0, 0.9*facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  //bump a non-adjacent pair of vertices in a random direction
  reload_mesh();
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}

 

