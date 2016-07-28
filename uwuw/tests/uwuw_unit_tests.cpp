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

  virtual void SetUp() {
    workflow_data = new UWUW(std::string(TEST_FILE));
  }
};

/*
 * Empty common setup function
 */
TEST_F(UWUWTest,SetUp)
{

}

/*
 * Test to make sure the total path is correct
 */
TEST_F(UWUWTest,filepath1)
{
  std::string filepath = "";
  EXPECT_NE(workflow_data->full_filepath,filepath);
  return;
}

/*
 * Test to make sure the total path is correct
 */
TEST_F(UWUWTest,filepath2)
{
  char current_path[FILENAME_MAX];
  getcwd(current_path,sizeof(current_path));
  std::string filepath(current_path);
  std::string filename(TEST_FILE);
  filepath+="/"+filename;
  EXPECT_EQ(workflow_data->full_filepath,filepath);
  return;
}

/*
 * Test of absolute path
 */
TEST_F(UWUWTest,filepath3)
{
  // get the full current path
  char current_path[FILENAME_MAX];
  // get the cwd
  getcwd(current_path,sizeof(current_path));
  std::string filepath(current_path);
  // full path to file
  filepath += "/"+std::string(TEST_FILE);
  // local uwuw class for this test only
  UWUW* wfd = new UWUW(filepath);
  EXPECT_EQ(wfd->full_filepath,filepath);
  return;
}

/*
 * Test of path with space in
 */
TEST_F(UWUWTest,filepath4)
{
  char current_path[FILENAME_MAX];
  // get the cwd
  getcwd(current_path,sizeof(current_path));
  // convert to std::string
  std::string filepath(current_path);
  std::string test_string = filepath+"/"+std::string(TEST_FILE)+" ";
  std::string correct_path = filepath+"/"+std::string(TEST_FILE);
  // local uwuw class for this test only
  UWUW* wfd = new UWUW(test_string);
  // expect filepath with last space removed
  EXPECT_EQ(wfd->full_filepath,correct_path);
  return;
}



/*
 * Test to make sure that the number of materials is correct
 */
TEST_F(UWUWTest,materiallibrary1)
{
  EXPECT_EQ(workflow_data->material_library.size(),2);
  return;
}

/*
 * Test to make sure that the materials were correctly loaded
 */
TEST_F(UWUWTest,materiallibrary2)
{
  // iterator for material library
  std::map<std::string,pyne::Material>::iterator it;
  it = workflow_data->material_library.begin();
  EXPECT_NE(it,workflow_data->material_library.end());
  return;
}

/*
 *  Test to make sure the the material can be read from any
 *  datapath in the file
 */
TEST_F(UWUWTest, material_datapath)
{
  // first we need to write some new materials
  pyne::comp_map nucvec;
  nucvec[pyne::nucname::id("H")] = 1.0;
  nucvec[pyne::nucname::id("Fe")] = 1.0;
  pyne::Material mat = pyne::Material(nucvec);
  mat.metadata["name"] = "Wet Steel";
  mat.write_hdf5("new_mat_test.h5","/materials"
                 ,"/nucid");

  pyne::Material mat2 = pyne::Material(nucvec);
  mat2.metadata["name"] = "Wet Steel 2";
  mat2.write_hdf5("new_mat_test.h5","/materials"
                  ,"/nucid");


  workflow_data->~UWUW();

  workflow_data = new UWUW(std::string("new_mat_test.h5"));

  // there should be 2 materials
  EXPECT_EQ(workflow_data->material_library.size(),2);
  return;
}

