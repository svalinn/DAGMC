//
// Patrick Shriwise
// June 2016
// This program is designed to run a set of tests on the make_watertight
// algorithm. This will be a stand-alone piece of code that uses MOAB to open,
// modify (break), and re-seal geometries input: cones.h5m output: pass/fail for
// each of the tests

#include "gtest/gtest.h"
#include "test_classes.hpp"

// Single Vertex Tests
TEST_F(MakeWatertightConeTest, SingleVertexMoveInXTest) {
  // bump in x direction
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.9 * facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}

TEST_F(MakeWatertightConeTest, SingleVertexMoveInYTest) {
  // bump in y direction
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.0, 0.9 * facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}

TEST_F(MakeWatertightConeTest, SingleVertexMoveInZTest) {
  // bump in z direction
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.0, 0.0, 0.9 * facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}

TEST_F(MakeWatertightConeTest, SingleVertexRandMoveTest) {
  // bump in random direction
  EXPECT_NO_THROW(result = rand_vert_bump(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
// Locked Pair Tests
TEST_F(MakeWatertightConeTest, LockedVertexPairMoveInXTest) {
  // bump locked pair in x direction
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.9 * facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, LockedVertexPairMoveInYTest) {
  // bump locked pair in y direction
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.0, 0.9 * facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, LockedVertexPairMoveInZTest) {
  // bump locked pair in z direction
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.0, 0.0, 0.9 * facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, LockedVertexPairRandMoveTest) {
  // bump locked pair in random diretion
  EXPECT_NO_THROW(result = locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
// Randomly Selected Vertex Pair Tests
TEST_F(MakeWatertightConeTest, RandomVertexPairMoveInXTest) {
  // bump random locked pair in x direction
  EXPECT_NO_THROW(result =
                      rand_locked_pair_bump(verts, 0.9 * facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, RandomVertexPairMoveInYTest) {
  // bump random locked pair in y direction
  EXPECT_NO_THROW(result =
                      rand_locked_pair_bump(verts, 0.0, 0.9 * facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, RandomVertexPairMoveInZTest) {
  // bump random locked pair in z direction
  EXPECT_NO_THROW(result =
                      rand_locked_pair_bump(verts, 0.0, 0.0, 0.9 * facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, RandomVertexPairRandMoveTest) {
  // bump random locked pair in random direction
  EXPECT_NO_THROW(result = rand_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
// Locked Pair With Vertex Between Tests
TEST_F(MakeWatertightConeTest, AdjacentPlusOneVertexPairMoveInXTest) {
  // bump pair of verices in the x direction
  EXPECT_NO_THROW(
      result = adjplone_locked_pair_bump(verts, 0.9 * facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, AdjacentPlusOneVertexPairMoveInYTest) {
  // bump pair of verices in the y direction
  EXPECT_NO_THROW(
      result = adjplone_locked_pair_bump(verts, 0.0, 0.9 * facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, AdjacentPlusOneVertexPairMoveInZTest) {
  // bump pair of verices in the z direction
  EXPECT_NO_THROW(
      result = adjplone_locked_pair_bump(verts, 0.0, 0.0, 0.9 * facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, AdjacentPlusOneVertexPairRandMoveTest) {
  // bump pair of verices in a random direction
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
// Randomly Selected Vertex Pair Tests
TEST_F(MakeWatertightConeTest, NonAdjacentVertexPairMoveInXTest) {
  // bump a non-adjacent pair of vertices in the x direction
  EXPECT_NO_THROW(
      result = nonadj_locked_pair_bump(verts, 0.9 * facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, NonAdjacentVertexPairMoveInYTest) {
  // bump a non-adjacent pair of vertices in the y direction
  EXPECT_NO_THROW(
      result = nonadj_locked_pair_bump(verts, 0.0, 0.9 * facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, NonAdjacentVertexPairMoveInZTest) {
  // bump a non-adjacent pair of vertices in the z direction
  EXPECT_NO_THROW(
      result = nonadj_locked_pair_bump(verts, 0.0, 0.0, 0.9 * facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightConeTest, NonAdjacentVertexPairRandMoveTest) {
  // bump a non-adjacent pair of vertices in a random direction
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
