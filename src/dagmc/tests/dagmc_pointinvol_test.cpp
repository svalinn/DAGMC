#include <gtest/gtest.h>

#include <iostream>

#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "DagMC.hpp"

using namespace moab;

using moab::DagMC;

std::shared_ptr<moab::DagMC> DAG;

static const char input_file[] = "test_geom.h5m";

class DagmcPointInVolTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    // Create new DAGMC instance
    DAG = std::make_shared<moab::DagMC>();
    // Load mesh from file
    rloadval = DAG->load_file(input_file);
    assert(rloadval == moab::MB_SUCCESS);
    // Create the OBB
    rval = DAG->init_OBBTree();
    assert(rval == moab::MB_SUCCESS);
  }
  virtual void TearDown() {}
 protected:
  moab::ErrorCode rloadval;
  moab::ErrorCode rval;
};

TEST_F(DagmcPointInVolTest, dagmc_setup_test) {
  ErrorCode rval = DAG->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = DAG->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcPointInVolTest, dagmc_point_in) {
  int result = 0;
  int expected_result = 1;
  double xyz[3] = {0.0, 0.0, 0.0};
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  ErrorCode rval = DAG->point_in_volume(vol_h, xyz, result);
  EXPECT_EQ(rval, MB_SUCCESS);
  EXPECT_EQ(expected_result, result);
}

int dagmc_point_in_vol_dir(double origin[3], double dir[3], int vol_idx) {
  int result = 0;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double xyz[3];
  double next_surf_dist;
  EntityHandle next_surf;

  // normalise the direction vector
  double dir_norm = (dir[0] * dir[0]) + (dir[1] * dir[1]) + (dir[2] * dir[2]);
  dir[0] = dir[0] / sqrt(dir_norm);
  dir[1] = dir[1] / sqrt(dir_norm);
  dir[2] = dir[2] / sqrt(dir_norm);

  ErrorCode rval = DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist);
  EXPECT_EQ(rval, MB_SUCCESS);

  xyz[0] = origin[0] + (next_surf_dist * dir[0]);
  xyz[1] = origin[1] + (next_surf_dist * dir[1]);
  xyz[2] = origin[2] + (next_surf_dist * dir[2]);

  std::cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;

  rval = DAG->point_in_volume(vol_h, xyz, result, dir);
  EXPECT_EQ(rval, MB_SUCCESS);
  return result;
}

TEST_F(DagmcPointInVolTest, dagmc_point_in_vol_1) {
  double dir[3] = {-1.0, 0.0, 0.0};
  double origin[3] = {0.0, 0.0, 0.0};
  int vol_idx = 1;
  int expected_result = 1;

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);
  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_in_vol_2) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {1.0, 0.0, 0.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_in_vol_3) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {0.0, -1.0, 0.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_in_vol_4) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {0.0, 1.0, 0.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_in_vol_5) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {0.0, 0.0, -1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_in_vol_6) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {0.0, 0.0, 1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_on_corner_1) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {1.0, 1.0, 1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_on_corner_2) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {-1.0, 1.0, 1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_on_corner_3) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {1.0, 1.0, -1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_on_corner_4) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {-1.0, 1.0, -1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_on_corner_5) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {1.0, -1.0, 1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_on_corner_6) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {-1.0, -1.0, 1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_on_corner_7) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {1.0, -1.0, -1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}

TEST_F(DagmcPointInVolTest, dagmc_point_on_corner_8) {
  int expected_result = 1;
  int vol_idx = 1;
  double dir[3] = {-1.0, -1.0, -1.0};
  double origin[3] = {0.0, 0.0, 0.0};

  int result = dagmc_point_in_vol_dir(origin, dir, vol_idx);

  EXPECT_EQ(expected_result, result);
}
