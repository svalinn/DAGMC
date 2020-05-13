#include <iostream>

#include <gtest/gtest.h>

#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "moab/GeomQueryTool.hpp"
#include "DagMC.hpp"

using namespace moab;

using moab::DagMC;

std::shared_ptr<moab::DagMC> DAG;

static const char input_file[] = "test_geom.h5m";
double eps = 1.0e-6;

class DagmcRayFireTest : public ::testing::Test {
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

TEST_F(DagmcRayFireTest, dagmc_setup_test) {
  ErrorCode rval = DAG->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = DAG->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcRayFireTest, dagmc_origin_face_rayfire) {
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {-1.0, 0.0, 0.0};
  double origin[3] = {0.0, 0.0, 0.0};
  double next_surf_dist;
  EntityHandle next_surf;
  DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist);
  double expected_next_surf_dist = 5.0;
  EXPECT_NEAR(expected_next_surf_dist, next_surf_dist, eps);
}

TEST_F(DagmcRayFireTest, dagmc_outside_face_rayfire) {
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0};  // ray along x direction
  double origin[3] = {-10.0, 0.0, 0.0};  // origin at -10 0 0
  double next_surf_dist;
  EntityHandle next_surf;
  DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist);
  std::cout << next_surf_dist << std::endl;
  double expected_next_surf_dist = 15.0;
  EXPECT_NEAR(expected_next_surf_dist, next_surf_dist, eps);
}

TEST_F(DagmcRayFireTest, dagmc_outside_face_rayfire_orient_exit) {
  DagMC::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0};  // ray along x direction
  double origin[3] = {-10.0, 0.0, 0.0};  // origin at -10 0 0
  double next_surf_dist;
  EntityHandle next_surf;
  DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
  std::cout << next_surf_dist << std::endl;
  double expected_next_surf_dist = 15.0;
  EXPECT_NEAR(expected_next_surf_dist, next_surf_dist, eps);
}

TEST_F(DagmcRayFireTest, dagmc_outside_face_rayfire_orient_entrance) {
  DagMC::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0};  // ray along x direction
  double origin[3] = {-10.0, 0.0, 0.0};  // origin at -10 0 0
  double next_surf_dist;
  EntityHandle next_surf;
  DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0.0,
                -1);
  std::cout << next_surf_dist << std::endl;
  double expected_next_surf_dist = 5.0;
  EXPECT_NEAR(expected_next_surf_dist, next_surf_dist, eps);
}

TEST_F(DagmcRayFireTest, dagmc_outside_face_rayfire_history_fail) {
  DagMC::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0};  // ray along x direction
  double origin[3] = {-10.0, 0.0, 0.0};  // origin at -10 0 0
  double xyz[3];
  double next_surf_dist;
  EntityHandle next_surf;

  history.reset();

  // ray fired exactly along boundary shared by 2 facets on a single surface,
  // needs two ray_fires to cross, this is expected and ok

  // first ray fire with history
  DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
  // second ray fire with history
  DAG->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist, &history, 0, 1);
  // this fire should hit graveyard, i.e. next_surf = 0
  DAG->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist, &history, 0, 1);

  // using history with this geom, there should be no next surface, i.e. 0
  EntityHandle ZERO = 0;
  EXPECT_EQ(ZERO, next_surf);
}

TEST_F(DagmcRayFireTest, dagmc_outside_face_rayfire_history) {
  DagMC::RayHistory history;
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  double dir[3] = {1.0, 0.0, 0.0};  // ray along x direction
  double origin[3] = {-10.0, 0.0, 0.0};  // origin at -10 0 0
  double xyz[3];
  double next_surf_dist;
  EntityHandle next_surf;

  history.reset();
  // first ray fire with history
  DAG->ray_fire(vol_h, origin, dir, next_surf, next_surf_dist, &history, 0, 1);
  std::cout << next_surf << " " << history.size() << std::endl;
  // second ray fire with history

  xyz[0] = origin[0] + (next_surf_dist * dir[0]);
  xyz[1] = origin[1] + (next_surf_dist * dir[1]);
  xyz[2] = origin[2] + (next_surf_dist * dir[2]);

  // ray fired execacyl

  DAG->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist, &history, 0, 1);

  DAG->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist, &history, 0, 1);

  // using history with this geom, there should be no next surface, i.e. 0
  EntityHandle ZERO = 0;
  EXPECT_EQ(ZERO, next_surf);
}
