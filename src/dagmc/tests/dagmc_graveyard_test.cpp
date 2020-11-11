#include <gtest/gtest.h>

#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "DagMC.hpp"

#include <iostream>

using namespace moab;

using moab::DagMC;

moab::DagMC* DAG;

static const char input_file[] = "test_dagmc.h5m";

class DagmcGraveyardTest : public ::testing::Test {
 protected:
  virtual void SetUp() {}
  virtual void TearDown() {}
};


class GraveyardTest : public::testing::Test {
  protected:

  virtual void SetUp() override {}
  virtual void TearDown() {}
};

TEST_F(DagmcGraveyardTest, dagmc_load_file) {
  DAG = new DagMC();
  ErrorCode rval = DAG->load_file(input_file); // open the Dag file
  EXPECT_EQ(rval, MB_SUCCESS);

  rval = DAG->init_OBBTree();
  EXPECT_EQ(rval, MB_SUCCESS);

  EXPECT_FALSE(DAG->has_graveyard());

  rval = DAG->create_graveyard();
  EXPECT_EQ(rval, MB_SUCCESS);

  EXPECT_TRUE(DAG->has_graveyard());
}
