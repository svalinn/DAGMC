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

#define TEST_FILE "test_uwuw.h5m"

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
 * Test to make sure the total path is correct
 */
TEST_F(UWUWTest,filepath1) {
  std::string filepath = "";
  EXPECT_NE(workflow_data->full_filepath,filepath);
  return;
}

/*
 * Test to make sure the total path is correct
 */
TEST_F(UWUWTest,filepath2) {
  char current_path[FILENAME_MAX];
  getcwd(current_path,sizeof(current_path));
  std::string filepath(current_path);
  std::string filename(TEST_FILE);
  filepath+="/"+filename;
  EXPECT_EQ(workflow_data->full_filepath,filepath);
  return;
}

/*
 * Test to make sure that the number of materials is correct
 */
TEST_F(UWUWTest,materiallibrary1) {
  EXPECT_EQ(workflow_data->material_library.size(),2);
  return;
}

/*
 * Test to make sure that the materials were correctly loaded
 */
TEST_F(UWUWTest,materiallibrary2) {
  // iterator for material library
  std::map<std::string,pyne::Material>::iterator it;
  it = workflow_data->material_library.begin();
  EXPECT_NE(it,workflow_data->material_library.end());
  return;
}



