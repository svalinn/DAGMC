#include "overlap_check_test.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

#include "moab/Core.hpp"
#include "overlap.hpp"

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
  virtual void SetFilename() override { filename = "overlap.h5m"; }
};

TEST_F(OverlappingVolumesTest, test1) {
  OverlapMap olaps;
  // triangle vertex locations only
  ErrorCode rval = check_instance_for_overlaps(MBI, olaps);
  EXPECT_EQ(rval, MB_SUCCESS);

  // we expect one overlap for this model between volmes 1 and 2
  EXPECT_EQ(olaps.size(), 1);
  std::set<int> expected_set = {1, 2};
  EXPECT_EQ(expected_set, olaps.begin()->first);
}

class NonOverlappingVolumesTest : public OverlapTest {
  virtual void SetFilename() override { filename = "no_overlap.h5m"; }
};

TEST_F(NonOverlappingVolumesTest, test2) {
  OverlapMap olaps;
  // triangle vertex locations only
  ErrorCode rval = check_instance_for_overlaps(MBI, olaps);
  EXPECT_EQ(rval, MB_SUCCESS);

  EXPECT_EQ(olaps.size(), 0);

  // test a few points along trianlge edges too
  rval = check_instance_for_overlaps(MBI, olaps, 2);
  EXPECT_EQ(rval, MB_SUCCESS);

  EXPECT_EQ(olaps.size(), 0);
}

class NonOverlappingImprintedVolumesTest : public OverlapTest {
  virtual void SetFilename() override { filename = "no_overlap_imp.h5m"; }
};

TEST_F(NonOverlappingImprintedVolumesTest, test3) {
  OverlapMap olaps;
  // triangle vertex locations only
  ErrorCode rval = check_instance_for_overlaps(MBI, olaps);
  EXPECT_EQ(rval, MB_SUCCESS);

  EXPECT_EQ(olaps.size(), 0);

  // test a few points along triangle edges too
  rval = check_instance_for_overlaps(MBI, olaps, 2);
  EXPECT_EQ(rval, MB_SUCCESS);

  EXPECT_EQ(olaps.size(), 0);
}

class EnclosedVolumeTest : public OverlapTest {
  virtual void SetFilename() override { filename = "enclosed.h5m"; }
};

TEST_F(EnclosedVolumeTest, test1) {
  OverlapMap olaps;
  // triangle vertex locations only
  ErrorCode rval = check_instance_for_overlaps(MBI, olaps);
  EXPECT_EQ(rval, MB_SUCCESS);

  // we expect one overlap for this model between volmes 1 and 2
  EXPECT_EQ(olaps.size(), 1);
  std::set<int> expected_set = {1, 2};
  EXPECT_EQ(expected_set, olaps.begin()->first);
}

class SmallOverlapTest : public OverlapTest {
  virtual void SetFilename() { filename = "small_overlap.h5m"; }
};

TEST_F(SmallOverlapTest, test1) {
  OverlapMap olaps;
  // triangle vertex locations only
  ErrorCode rval = check_instance_for_overlaps(MBI, olaps);
  EXPECT_EQ(rval, MB_SUCCESS);

  // we expect one overlap for this model between volmes 1 and 2
  EXPECT_EQ(olaps.size(), 1);
  std::set<int> expected_set = {1, 2};
  EXPECT_EQ(expected_set, olaps.begin()->first);
}
