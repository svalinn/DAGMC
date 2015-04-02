//  DagSolid_test.cpp
#include <gtest/gtest.h>
#include <cmath>
#include <cassert>

#include <iostream>
#include <unistd.h>
#include <stdio.h>

#include "uwuw.hpp"
#include "../pyne/pyne.h"

UWUW *workflow_data;

#define TEST_FILE "mat_lib.h5"

class UWUWTest : public ::testing::Test
{
  protected:

  virtual void SetUp()
  {
    workflow_data = new UWUW(TEST_FILE);
  }
};

/*
 * Empty common setup function
 */
TEST_F(UWUWTest,SetUp) {

}

/*
 * Test to make sure that the number of tallies is correct
 */
TEST_F(UWUWTest,tally_library_1) {
  EXPECT_EQ(workflow_data->tally_library.size(),0);
  return;
}
