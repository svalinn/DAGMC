//  DagSolid_test.cpp
#include <gtest/gtest.h>
#include "uwuw.hpp"
#include "../pyne/pyne.h"

#define TEST_FILE "mat_lib.h5"

namespace
{


class UWUWTest : public ::testing::Test
{
 protected:

  UWUWTest() {}
  virtual ~UWUWTest() {}

  virtual void SetUp() {
    workflow_data = new UWUW(std::string(TEST_FILE));
  }
  UWUW *workflow_data;
};

/*
 * Test to make sure that the number of tallies is correct
 */
TEST_F(UWUWTest,TallyLibraryEmpty)
{
  EXPECT_EQ(workflow_data->tally_library.size(),0);
  return;
}

TEST_F(UWUWTest,MaterialLibrarySomeMaterials)
{
  EXPECT_NE(workflow_data->material_library.size(),0);
  return;
}

TEST_F(UWUWTest,MaterialLibraryCorrectNumber)
{
  EXPECT_EQ(workflow_data->material_library.size(),12);
  return;
}

} // namespace
