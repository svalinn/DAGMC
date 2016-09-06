//
// Patrick Shriwise
// September 2016
// This program is designed to run a set of tests on the make_watertight algorithm.
// This will be a stand-alone piece of code that uses MOAB to open, modify
// (break), and re-seal geometries
// input: cyl.h5m file (found in ../make_watertight/test/)
// output: pass/fail for each of the tests

#include "test_classes.hpp"
#include "gtest/gtest.h"


//Make sure the sphere is not deleted when sealing
TEST_F(MakeWatertightNoCurveSphereTest, NoCurveSphereDeletionCheck)
{
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));

  //make sure the sphere wasn't deleted
  int entity_dimension = 3;
  int num_ents_expected = 1;
  EXPECT_NO_THROW(result = check_num_ents(entity_dimension, num_ents_expected));
  EXPECT_TRUE(result == moab::MB_SUCCESS);
}
