#include <gtest/gtest.h>

// std includes
#include <iostream>
#include <string>
// moab includes
#include "moab/ProgOptions.hpp"
#include "moab/ScdInterface.hpp"
#include "DagMC.hpp"


using namespace moab;


class DagmcPrecondTest : public ::testing::Test {

protected:

  DagMC* dagmc;
  SignedDistanceField* box;
  ErrorCode rval;
  std::string filename;

  void set_filename() {
    // sphere of radius 5, centered on the origin
    filename = "sphere_rad5.h5m";
  };
  
  virtual void SetUp() {
    
    set_filename();

    dagmc = new DagMC();
    
    //create a new dagmc instance and load the file
    rval = dagmc->load_file(filename.c_str());
    MB_CHK_SET_ERR_RET(rval,"Could not load the preconditioner test file");

    //get the ScdBox of the first (and only) volume
    Range vols;
    const int three = 3;
    const void* const three_val[] = {&three};
    Tag geomTag = dagmc->geom_tag();
    rval = dagmc->moab_instance()->get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag, three_val, 1, vols );
    MB_CHK_SET_ERR_RET(rval,"Could not get the model volumes");

    //create a new group
    EntityHandle sdf_group_set;
    rval = dagmc->moab_instance()->create_meshset(0, sdf_group_set);
    MB_CHK_SET_ERR_RET(rval, "Failed to create a new group for storing signed distance field volumes");

    Tag cat_tag;
    rval = dagmc->moab_instance()->tag_get_handle(CATEGORY_TAG_NAME, cat_tag);
    MB_CHK_SET_ERR_RET(rval, "Failed to get the category tag handle");

    std::string group = "Group";
    group.resize(32);
    rval = dagmc->moab_instance()->tag_set_data(cat_tag, &sdf_group_set, 1, group.c_str());
    MB_CHK_SET_ERR_RET(rval, "Failed to set the category tag of the signed distance field group set");
  
    std::string sdf = "sdf_0.5";
    rval = dagmc->moab_instance()->tag_set_data(dagmc->name_tag(), &sdf_group_set, 1, (void*)sdf.c_str());
    MB_CHK_SET_ERR_RET(rval, "Failed to set the name tag of the signed distance field group set");

    //there should only be one volume before we init the OBB tree
    if(1 != vols.size()) {
      MB_CHK_SET_ERR_RET(MB_FAILURE, "Incorrect number of volumes found after loading the model");
    }

    rval = dagmc->moab_instance()->add_entities(sdf_group_set, vols);
    MB_CHK_SET_ERR_RET(rval, "Failed to add the volume to the signed distance field group set");

    //initalize the OBBTree (and the preconditioning datastructs)
    rval = dagmc->init_OBBTree();
    MB_CHK_SET_ERR_RET(rval,"Could not initialize the OBBTree");

    rval = dagmc->build_preconditioner();
    MB_CHK_SET_ERR_RET(rval,"Could not construct preconditioner.");

  }

  virtual void TearDown() {
    delete dagmc;
    
  }
};

TEST_F(DagmcPrecondTest, precondition_ray_fire_test) {

  // get the volume (there should be only one, index starts at one)
  EntityHandle vol = dagmc->entity_by_index(3, 1);

  double ray_start[3] = {0.0, 0.0, 0.0};

  double ray_dir[3] = {1.0, 0.0, 0.0};

  // for a sphere of radius 5 and a SDF with size 0.1,
  // this should be a preconditionable length
  double ray_len = 1;

  EntityHandle next_surf;
  double next_surf_dist;
  bool preconditioned;

  rval = dagmc->precondition_ray_fire(vol, ray_start, ray_dir, ray_len, next_surf, next_surf_dist, preconditioned);
  MB_CHK_SET_ERR_RET(rval, "Failed to precondition ray");

  EXPECT_TRUE(preconditioned);

  // increase ray distance to an un-preconditionable length
  ray_len = 10.0;

  rval = dagmc->precondition_ray_fire(vol, ray_start, ray_dir, ray_len, next_surf, next_surf_dist, preconditioned);
  MB_CHK_SET_ERR_RET(rval, "Failed to precondition ray");

  EXPECT_TRUE(!preconditioned);

  // increase ray distance to an un-preconditionable length
  ray_len = 4.99;

  rval = dagmc->precondition_ray_fire(vol, ray_start, ray_dir, ray_len, next_surf, next_surf_dist, preconditioned);
  MB_CHK_SET_ERR_RET(rval, "Failed to precondition ray");

  EXPECT_TRUE(!preconditioned);

}
