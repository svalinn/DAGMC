//  DagSolid_test.cpp
#include <gtest/gtest.h>
#include "uwuw.hpp"
#include "../pyne/pyne.h"

#define TEST_FILE "mat_lib.h5"

namespace {


class UWUWTest : public ::testing::Test
{
  protected:

  UWUWTest() {}
  virtual ~UWUWTest(){}

  virtual void SetUp()
  {
    workflow_data = new UWUW("mat_lib.h5");
  }
  UWUW *workflow_data;
};

/*
 * Test to make sure that the number of tallies is correct
 */
TEST_F(UWUWTest,tallylibrary1) {
  EXPECT_EQ(workflow_data->tally_library.size(),0);
  return;
}
} // namespace
