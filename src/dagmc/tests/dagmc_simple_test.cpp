#include <gtest/gtest.h>

#include <iostream>

#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"

using namespace moab;

using moab::DagMC;

moab::DagMC* DAG;

static const char input_file[] = "test_geom.h5m";

class DagmcSimpleTest : public ::testing::Test {
 protected:
  virtual void SetUp() {}
  virtual void TearDown() {}
};

TEST_F(DagmcSimpleTest, dagmc_load_file) {
  DAG = new DagMC();
  ErrorCode rval = DAG->load_file(input_file);  // open the Dag file
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc) {
  /* 1 - Test with external moab, load file in DAGMC*/
  // make new moab core
  std::shared_ptr<Interface> mbi = std::make_shared<Core>();
  // make new dagmc into that moab
  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>(mbi);

  ErrorCode rval;

  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc_via_moab) {
  /* 2 - Test with external moab, load file in MOAB*/
  // load the file into moab rather than dagmc
  ErrorCode rval;

  std::shared_ptr<Interface> mbi = std::make_shared<Core>();
  rval = mbi->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>(mbi);
  rval = dagmc->load_existing_contents();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc_internal) {
  /* 3 - Test with internal moab, load file in DAG*/
  // make new dagmc into that moab
  ErrorCode rval;

  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>();
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc_build_obb) {
  /* 1 - Test with external moab, load file in DAGMC*/
  // make new moab core
  ErrorCode rval;

  std::shared_ptr<Interface> mbi = std::make_shared<Core>();
  // make new dagmc into that moab
  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>(mbi);

  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc_via_moab_build_obb) {
  /* 2 - Test with external moab, load file in MOAB*/
  // load the file into moab rather than dagmc
  ErrorCode rval;

  std::shared_ptr<Interface> mbi = std::make_shared<Core>();
  rval = mbi->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>(mbi);
  rval = dagmc->load_existing_contents();
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_load_file_dagmc_internal_build_obb) {
  /* 3 - Test with internal moab, load file in DAG*/
  // make new dagmc into that moab
  ErrorCode rval;

  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>();
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_test_obb_retreval) {
  // make new dagmc
  std::cout << "test_obb_retreval" << std::endl;

  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>();

  ErrorCode rval;
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);

  // write the file
  rval = dagmc->write_mesh("fcad", 4);

  dagmc.reset(new DagMC());
  rval = dagmc->load_file("fcad");
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);

  // delete the fcad file
  remove("fcad");
}

TEST_F(DagmcSimpleTest, dagmc_build_obb) {
  ErrorCode rval = DAG->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);
}

TEST_F(DagmcSimpleTest, dagmc_num_vols) {
  int expect_num_vols = 2;
  int num_vols = DAG->num_entities(3);
  EXPECT_EQ(expect_num_vols, num_vols);
}

TEST_F(DagmcSimpleTest, dagmc_point_in) {
  int result = 0;
  int expect_result = 1;
  int vol_idx = 1;
  double xyz[3] = {0.0, 0.0, 0.0};
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  ErrorCode rval = DAG->point_in_volume(vol_h, xyz, result);
  EXPECT_EQ(rval, MB_SUCCESS);
  EXPECT_EQ(expect_result, result);
}

TEST_F(DagmcSimpleTest, dagmc_test_obb_retreval_rayfire) {
  // make new dagmc
  std::cout << "test_obb_retreval and ray_fire" << std::endl;

  std::shared_ptr<DagMC> dagmc = std::make_shared<DagMC>();

  ErrorCode rval;
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);

  // write the file
  rval = dagmc->write_mesh("fcad", 4);

  // now create new DAGMC
  dagmc.reset(new DagMC());
  rval = dagmc->load_file("fcad");
  EXPECT_EQ(rval, MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);

  // delete the fcad file
  remove("fcad");

  // now perform full ray fire
  double eps = 1.e-6;
  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {0.0, 0.0, 0.0};
  double dir[3] = {0.0, 0.0, 1.0};
  EntityHandle next_surf;
  double next_surf_dist;
  double expect_next_surf_dist = 5.0;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);

  rval = DAG->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist);
  EXPECT_EQ(rval, MB_SUCCESS);
  EXPECT_NEAR(expect_next_surf_dist, next_surf_dist, eps);
}

TEST_F(DagmcSimpleTest, dagmc_rayfire) {
  const double eps = 1e-6;  // epsilon for test, faceting tol?

  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {0.0, 0.0, 0.0};
  double dir[3] = {0.0, 0.0, 1.0};
  EntityHandle next_surf;
  double next_surf_dist;
  double expect_next_surf_dist = 5.0;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);

  ErrorCode rval = DAG->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist);
  EXPECT_EQ(rval, MB_SUCCESS);
  EXPECT_NEAR(expect_next_surf_dist, next_surf_dist, eps);
}

TEST_F(DagmcSimpleTest, dagmc_closest_to) {
  const double eps = 1e-6;  // epsilon for test, faceting tolerance

  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {-6.0, 0.0, 0.0};
  double distance;  // distance from point to nearest surface
  double expect_distance = 1.0;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);

  ErrorCode rval = DAG->closest_to_location(vol_h, xyz, distance);
  EXPECT_EQ(rval, MB_SUCCESS);
  // distance should be 1.0 cm
  EXPECT_NEAR(expect_distance, distance, eps);
}

TEST_F(DagmcSimpleTest, dagmc_test_boundary) {
  int vol_idx = 1;
  EntityHandle vol_h = DAG->entity_by_index(3, vol_idx);
  int surf_idx = 1;
  EntityHandle surf_h = DAG->entity_by_index(2, surf_idx);

  double xyz[3] = {0.0, 0.0, 5.0};
  double dir[3] = {0.0, 0.0, 1.0};
  int result;
  int expect_result = 0;

  ErrorCode rval = DAG->test_volume_boundary(vol_h, surf_h, xyz, dir, result);
  EXPECT_EQ(rval, MB_SUCCESS);
  // check ray leaving volume
  EXPECT_EQ(expect_result, result);
}
