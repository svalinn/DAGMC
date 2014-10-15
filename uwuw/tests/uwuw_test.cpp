//  DagSolid_test.cpp
#include <gtest/gtest.h>
#include <cmath>
#include <cassert>

#include <iostream>
#include "uwuw.hpp"

class UWUWTest : public ::testing::Test
{
  protected:

  virtual void SetUp()
  {
    UWUW workflow_data = UWUW("test.h5m");
  }

  protected:
  
  UWUW workflow_data;

};

/*
 * Empty common setup function
 */
TEST_F(UWUWTest,SetUp) {

}

/*
 * Test to make sure the total path is correct
 */
TEST_F(UWUWTest,filepath) {
  std::string filepath = "test";
  EXPECT_EQ(workflow_data.full_filepath,filepath);
  return;
}

