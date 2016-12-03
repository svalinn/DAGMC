// FluDAG/src/test/test_FlukaFuncs.cpp

#include <gtest/gtest.h>

#include "DagMC.hpp"
#include "moab/Interface.hpp"
#include "fluka_funcs.h"


#include <cmath>
#include <cassert>

moab::DagMC *DAG = new moab::DagMC();

int num_slab_vols = 12;

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class FluDAGTest : public ::testing::Test
{
 protected:

  // initalize variables for each test
  virtual void SetUp() {
    // Default h5m file for testing
    std::string infile = "slabs.h5m";

    rloadval = DAG->load_file(infile.c_str());
    assert(rloadval == moab::MB_SUCCESS);

    // DAG call to initialize geometry
    rval = DAG->init_OBBTree();
    assert (rval == moab::MB_SUCCESS);

    // Initialize point and dir
    point[0] = 0.0;
    point[1] = 0.0;
    point[2] = 0.0;

    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 0.0;

    oldReg = 1;
    // How far the particle would be expected to go if no physics boundaries
    propStep = 0.0;
    safety   = 0.0;

    // Direction cosine for component headed from center of cube to a corner
    dir_norm = 1.0/sqrt(3);
  }

 protected:

  moab::ErrorCode rloadval;
  moab::ErrorCode rval;

  // Position
  double point[3];

  // Direction Vector
  double dir[3];

  int      oldReg;
  double propStep;
  double safety;
  double retStep;
  int    newReg;

  double dir_norm;
};

//---------------------------------------------------------------------------//
// void g_fire(int& oldRegion, double point[], double dir[],
//              double &propStep, double& retStep, double& safety, int& newRegion)
//---------------------------------------------------------------------------//
// oldRegion - the region of the particle's current coordinates
// point     - the particle's coordinate location vector
// dir       - the direction vector of the particle's current path (ray)
// propStep  - ??
// retStep   - returned as the distance from the particle's current location, along its ray, to the next boundary
// newRegion - gotten from the value returned by DAG->next_vol
// newRegion is gotten from the volue returned by DAG->next_vol
// void g_fire(int& oldRegion, double point[], double dir[], double &propStep, double& retStep, double& safety, int& newRegion)

//---------------------------------------------------------------------------//
// Test setup outcomes
TEST_F(FluDAGTest, SetUp)
{
  EXPECT_EQ(moab::MB_SUCCESS, rloadval);

  // DAG call to initialize geometry
  EXPECT_EQ(moab::MB_SUCCESS, rval);

  int num_vols = DAG->num_entities(3);
  std::cout << "Number of regions is " << num_vols << std::endl;
  EXPECT_EQ(num_slab_vols, num_vols);

  std::vector< std::string > keywords;
  rval = DAG->detect_available_props( keywords );
  EXPECT_EQ(moab::MB_SUCCESS, rval);
  rval = DAG->parse_properties( keywords );
  EXPECT_EQ(moab::MB_SUCCESS, rval);

  int ret, volume;
  for (unsigned i=1; i<=num_slab_vols; i++) {
    int id = DAG->id_by_index(3, i);
    std::cout << "Vol " << i << ", id = " << id << std::endl;

    moab::EntityHandle eh = DAG->entity_by_index(3,i);
    rval = DAG->point_in_volume(eh, point, ret);
    EXPECT_EQ(moab::MB_SUCCESS, rval);
    if (ret == 1) {
      volume = i;
      std::cout << "\tPoint is in this volume!" << std::endl;
    }
  }
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: WrapperTest
//---------------------------------------------------------------------------//
// Test default propStep value of 0.0 is not right
TEST_F(FluDAGTest, GFireBadPropStep)
{
  std::cout << "Calling g_fire. Start in middle of leftmost 10x10x10 cube" << std::endl;
  oldReg   = 2;
  point[2] = 5.0;
  dir[2]   = 1.0;

  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_EQ(0.0, retStep);
}
//---------------------------------------------------------------------------//
// Test distance to next surface or vertex for various directions
TEST_F(FluDAGTest, GFireGoodPropStep)
{
  // std::cout << "Calling g_fire. Start in middle of leftmost cube" << std::endl;
  oldReg   = 2;
  point[2] = 5.0;
  // Set prepStep to something more realistic than 0.0
  propStep = 1e38;

  // +z direction
  dir[2]   = 1.0;
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  std::cout << "newReg is " << newReg << std::endl;
  // Start in middle of 10x10x10 cube, expect 10/2 to be dist. to next surface
  EXPECT_EQ(5.0, retStep);
  // -z direction
  dir[2] = -dir[2];
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_EQ(5.0, retStep);
  // +y direction
  dir[2] = 0.0;
  dir[1] = 1.0;
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_EQ(5.0, retStep);
  // -y direction
  dir[1] = -dir[1];
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_EQ(5.0, retStep);
  // +x direction
  dir[1] = 0.0;
  dir[0] = 1.0;
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_EQ(5.0, retStep);
  // -x direction
  dir[0] = -dir[0];
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_EQ(5.0, retStep);

  // +++
  dir[0] = +dir_norm;
  dir[1] = +dir_norm;
  dir[2] = +dir_norm;
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_NEAR(8.660254, retStep, 1e-6);

  // ++-
  // Not Lost Particle!
  dir[0] = +dir_norm;
  dir[1] = +dir_norm;
  dir[2] = -dir_norm;
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_NEAR(8.660254, retStep, 1e-6);

  // +-+
  dir[0] = +dir_norm;
  dir[1] = -dir_norm;
  dir[2] = +dir_norm;
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_NEAR(8.660254, retStep, 1e-6);

  // +--
  // Not Lost Particle!
  dir[0] = +dir_norm;
  dir[1] = -dir_norm;
  dir[2] = -dir_norm;
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_NEAR(8.660254, retStep, 1e-6);

  // -++
  dir[0] = -dir_norm;
  dir[1] = +dir_norm;
  dir[2] = +dir_norm;
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_DOUBLE_EQ(5.0/dir_norm, retStep);

  // -+-
  // Not Lost Particle!
  dir[0] = -dir_norm;
  dir[1] = +dir_norm;
  dir[2] = -dir_norm;
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_DOUBLE_EQ(5.0/dir_norm, retStep);

  // --+
  dir[0] = -dir_norm;
  dir[1] = -dir_norm;
  dir[2] = +dir_norm;
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_DOUBLE_EQ(5.0/dir_norm, retStep);

  // ---
  // Not Lost Particle!
  dir[0] = -dir_norm;
  dir[1] = -dir_norm;
  dir[2] = -dir_norm;
  g_fire(oldReg, point, dir, propStep, retStep, safety, newReg);
  EXPECT_DOUBLE_EQ(5.0/dir_norm, retStep);
}

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class FluDAGMetaDataTest : public ::testing::Test
{
};


//---------------------------------------------------------------------------//
// Test setup outcomes
TEST_F(FluDAGMetaDataTest, CheckAssignMatsLegacy)
{
  delete DAG;
  DAG = new moab::DagMC();
  std::string infile = "test_geom_legacy.h5m";
  moab::ErrorCode error;
  error = DAG->load_file(infile.c_str()); // load the dag file takeing the faceting from h5m
  error = DAG->setup_impl_compl();
  error = DAG->setup_indices();
  fludag_write(infile,"lcad");

// expected values from the lcad file // only the cells
  const char* expected[] = {"*...+....1....+....2....+....3....+....4....+....5....+....6....+....7...",
                            "ASSIGNMA           9        1.",
                            "ASSIGNMA           9        2.",
                            "ASSIGNMA           9        3.",
                            "ASSIGNMA           9        4.",
                            "ASSIGNMA           3        5.",
                            "ASSIGNMA           3        6.",
                            "ASSIGNMA           3        7.",
                            "ASSIGNMA           3        8.",
                            "ASSIGNMA           1        9.",
                            "ASSIGNMA    BLCKHOLE       10.",
                            "ASSIGNMA      VACUUM       11."
                           };
  std::vector<std::string> expected_lcad(expected,expected+12);

  // now read the lcad file
  std::ifstream input;
  input.open("lcad");
  std::string line;
  std::vector<std::string> input_deck;
  while(!input.eof()) {
    std::getline(input,line);
    input_deck.push_back(line);
  }
  input.close();

  // for each line make sure the same
  for ( int i = 0 ; i < 11 ; i++ ) {
    EXPECT_EQ(expected_lcad[i],input_deck[i]);
  }
  // delete the lcad file
  std::remove("lcad");

}