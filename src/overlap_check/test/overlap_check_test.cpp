//
// Patrick Shriwise
// Sept 2019
// This program tests for overlaps in some known test models - some containing
// overlapping volume, some not.
// input: overlap.h5m
// output: pass/fail for each of the tests

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <algorithm>
#include <memory>
#include "moab/Core.hpp"

#include "overlap.hpp"

#include "overlap_check_test.hpp"

void OverlapTest::SetUp() {
  SetFilename();

  MBI = std::shared_ptr<Interface>(new Core());

  ErrorCode rval = MBI->load_file(filename.c_str());
  EXPECT_EQ(rval, MB_SUCCESS);
};

void OverlapTest::TearDown() {
  ErrorCode rval = MBI->delete_mesh();
  EXPECT_EQ(rval, MB_SUCCESS);
};

class OverlappingVolumesTest : public OverlapTest {
  virtual void SetFilename() { filename = "overlap.h5m"; }
};

TEST_F(OverlappingVolumesTest, test1) {
  OverlapMap olaps;
  ErrorCode rval = check_file_for_overlaps(MBI, olaps);
  EXPECT_EQ(rval, MB_SUCCESS);

  // we expect one overlap for this model between volmes 1 and 2
  EXPECT_EQ(olaps.size(), 1);
  std::set<int> expected_set = {1,2};
  EXPECT_EQ(expected_set, olaps.begin()->first);
}

class NonOverlappingVolumesTest : public OverlapTest {
  virtual void SetFilename() { filename = "no_overlap.h5m"; }
};

TEST_F(NonOverlappingVolumesTest, test2) {
  OverlapMap olaps;
  // tri verts only test
  ErrorCode rval = check_file_for_overlaps(MBI, olaps);
  EXPECT_EQ(rval, MB_SUCCESS);

  EXPECT_EQ(olaps.size(), 0);

  // test tri edges too
  rval = check_file_for_overlaps(MBI, olaps, 2);
  EXPECT_EQ(rval, MB_SUCCESS);

  EXPECT_EQ(olaps.size(), 0);
}

class NonOverlappingImprintedVolumesTest : public OverlapTest {
  virtual void SetFilename() { filename = "no_overlap_imp.h5m"; }
};

TEST_F(NonOverlappingImprintedVolumesTest, test3) {
  OverlapMap olaps;
  // tri verts only test
  ErrorCode rval = check_file_for_overlaps(MBI, olaps);
  EXPECT_EQ(rval, MB_SUCCESS);

  EXPECT_EQ(olaps.size(), 0);

  // test tri edges too
  rval = check_file_for_overlaps(MBI, olaps, 2);
  EXPECT_EQ(rval, MB_SUCCESS);

  EXPECT_EQ(olaps.size(), 0);

}
