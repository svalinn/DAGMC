//
// Patrick Shriwise
// Sept 2019
// This program tests for overlaps in some known test models - some containing
// overlapping volume, some not.
// input: overlap.h5m
// output: pass/fail for each of the tests


#include "gtest/gtest.h"


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <algorithm>
#include <memory>
#include "moab/Core.hpp"

using namespace moab;
class OverlapCheckTest : public::testing::Test {
protected:
  void SetUp() {
    setFilename();

    MBI = (new Core());

    ErrorCode rval = MBI->load_file(filename.c_str());
    EXPECT_EQ(rval, MB_SUCCESS);
  };

  virtual void TearDown() {
    ErrorCode rval = MBI->delete_mesh();
    EXPECT_EQ(rval, MB_SUCCESS);
  };

  virtual void setFilename();

  std::string filename;

public:
  std::shared_ptr<Interface> MBI;
};

class OverlapFileTest : public OverlapCheckTest {

  virtual void setFilename() {
    filename = "overlap.h5m";
  }

};

TEST_F(OverlapCheckTest, OverlapCheckTest) {
  OverlapMap olaps;
  ErrorCode rval = check_file_for_overlaps(MBI, olaps);
  EXPECT_EQ(rval, MB_SUCCESS);

  EXPECT_EQ(olaps.size(), 0);
}
