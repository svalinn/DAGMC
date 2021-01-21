//
// Patrick Shriwise
// September 2013
// This program is designed to run a set of tests on the make_watertight
// algorithm. This will be a stand-alone piece of code that uses MOAB to open,
// modify (break), and re-seal geometries input: cyl.h5m file (found in
// ../make_watertight/test/) output: pass/fail for each of the tests

#include "gtest/gtest.h"
#include "test_classes.hpp"

// Single Vertex Tests
TEST_F(MakeWatertightCylinderTest, SingleVertexMoveInXTest) {
  // bump vert in x direction
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.9 * facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, SingleVertexMoveInYTest) {
  // bump vert in y direction
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.0, 0.9 * facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, SingleVertexMoveInZTest) {
  // bump vert in z direction
  EXPECT_NO_THROW(result = single_vert_bump(verts, 0.0, 0.0, 0.9 * facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, SingleVertexMoveInRTest) {
  // bump vert away from the center of the cylinder
  EXPECT_NO_THROW(result = single_vert_bump_R(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, SingleVertexMoveInThetaTest) {
  // bump vert along theta of the cylinder
  EXPECT_NO_THROW(result = theta_vert_bump(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, SingleVertexRandMoveTest) {
  // bump vert in a random direction
  EXPECT_NO_THROW(result = rand_vert_bump(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
// Locked Pair Tests
TEST_F(MakeWatertightCylinderTest, LockedVertexPairMoveInXTest) {
  // bump a pair of verts in the same amount in the x direction
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.9 * facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, LockedVertexPairMoveInYTest) {
  // bump a pair of verts in the same amount in the y direction
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.0, 0.9 * facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, LockedVertexPairMoveInZTest) {
  // bump a pair of verts in the same amount in the z direction
  EXPECT_NO_THROW(result = locked_pair_bump(verts, 0.0, 0.0, 0.9 * facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, LockedVertexPairMoveInRTest) {
  // bump a pair of verts in the same amount in the radial direction
  EXPECT_NO_THROW(result = locked_pair_bump_R(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, LockedVertexPairMoveInThetaTest) {
  // bump a pair of verts in the same amount in the theta direction
  EXPECT_NO_THROW(result = locked_pair_bump_theta(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, LockedVertexPairRandMoveTest) {
  // bump a pair of verts in the same amount in a random direction
  EXPECT_NO_THROW(result = locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
// Randomly Selected Vertex Pair Tests
TEST_F(MakeWatertightCylinderTest, RandomVertexPairMoveInXTests) {
  // move a randomly selected locked pair of vertices in the x direction
  EXPECT_NO_THROW(result =
                      rand_locked_pair_bump(verts, 0.9 * facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, RandomVertexPairMoveInYTest) {
  // move a randomly selected locked pair of vertices in the y direction
  EXPECT_NO_THROW(result =
                      rand_locked_pair_bump(verts, 0.0, 0.9 * facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, RandomVertexPairMoveInZTest) {
  // move a randomly selected locked pair of vertices in the z direction
  EXPECT_NO_THROW(result =
                      rand_locked_pair_bump(verts, 0.0, 0.0, 0.9 * facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, RandomVertexPairMoveInRTest) {
  // move a randomly selected locked pair of vertices in the radial direction
  EXPECT_NO_THROW(result = rand_locked_pair_bump_R(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, RandomVertexPairThetaMoveTest) {
  // move a randomly selected locked pair of vertices along the theta direction
  EXPECT_NO_THROW(result = rand_locked_pair_bump_theta(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, RandomVertexPairRandMoveTest) {
  // move a randomly selected locked pair of vertices in a random direction
  EXPECT_NO_THROW(result = rand_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
// Locked Pair With Vertex Between Tests
TEST_F(MakeWatertightCylinderTest, AdjacentPlusOneVertexPairMoveInXTest) {
  // bump a vertex pair with a single, stationary vertex between them in the x
  // direction
  EXPECT_NO_THROW(
      result = adjplone_locked_pair_bump(verts, 0.9 * facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, AdjacentPlusOneVertexPairMoveInYTest) {
  // bump a vertex pair with a single, stationary vertex between them in the y
  // direction
  EXPECT_NO_THROW(
      result = adjplone_locked_pair_bump(verts, 0.0, 0.9 * facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, AdjacentPlusOneVertexPairMoveInZTest) {
  // bump a vertex pair with a single, stationary vertex between them in the z
  // direction
  EXPECT_NO_THROW(
      result = adjplone_locked_pair_bump(verts, 0.0, 0.0, 0.9 * facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, AdjacentPlusOneVertexPairMoveInRTest) {
  // bump a vertex pair with a single, stationary vertex between them in the
  // radial direction
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump_R(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, AdjacentPlusOneVertexPairMoveInThetaTest) {
  // bump a vertex pair with a single, stationary vertex between them along the
  // theta direction
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump_theta(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, AdjacentPlusOneVertexPairRandMoveTest) {
  // bump a vertex pair with a single, stationary vertex between them in a
  // random direction
  EXPECT_NO_THROW(result = adjplone_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
// Randomly Selected Vertex Pair Tests
TEST_F(MakeWatertightCylinderTest, NonAdjacentVertexPairMoveInXTest) {
  // move two randomly selected, non-adjacent vertices in the x direction
  EXPECT_NO_THROW(
      result = nonadj_locked_pair_bump(verts, 0.9 * facet_tol, 0.0, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, NonAdjacentVertexPairMoveInYTest) {
  // move two randomly selected, non-adjacent vertices in the y direction
  EXPECT_NO_THROW(
      result = nonadj_locked_pair_bump(verts, 0.0, 0.9 * facet_tol, 0.0));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, NonAdjacentVertexPairMoveInZTest) {
  // move two randomly selected, non-adjacent vertices in the z direction
  EXPECT_NO_THROW(
      result = nonadj_locked_pair_bump(verts, 0.0, 0.0, 0.9 * facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, NonAdjacentVertexPairMoveInRTest) {
  // move two randomly selected, non-adjacent vertices in the radial direction
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump_R(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, NonAdjacentVertexPairMoveInThetaTest) {
  // move two randomly selected, non-adjacent vertices along the theta direction
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump_theta(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
TEST_F(MakeWatertightCylinderTest, NonAdjacentVertexPairRandMoveTest) {
  // move two randomly selected, non-adjacent vertices in a random direction
  EXPECT_NO_THROW(result = nonadj_locked_pair_bump_rand(verts, facet_tol));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
