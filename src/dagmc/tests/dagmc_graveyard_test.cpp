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
  EXPECT_EQ(MB_SUCCESS, rval);

  rval = DAG->init_OBBTree();
  EXPECT_EQ(MB_SUCCESS, rval);

  Range starting_vertices;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBVERTEX, starting_vertices);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range starting_triangles;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBTRI, starting_triangles);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range starting_sets;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBENTITYSET, starting_sets);
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_FALSE(DAG->has_graveyard());

  int n_vols = DAG->num_entities(3);

  rval = DAG->create_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_EQ(n_vols + 1, DAG->num_entities(3));

  EXPECT_TRUE(DAG->has_graveyard());

  rval = DAG->remove_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_FALSE(DAG->has_graveyard());

  rval = DAG->create_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_TRUE(DAG->has_graveyard());

  rval = DAG->remove_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_FALSE(DAG->has_graveyard());

  // checks to make sure we didn't accumulate any new data on the mesh
  Range ending_vertices;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBVERTEX, ending_vertices);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range ending_triangles;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBTRI, ending_triangles);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range ending_sets;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBENTITYSET, ending_sets);
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_EQ(starting_vertices.size(), ending_vertices.size());
  EXPECT_EQ(0, subtract(starting_vertices, ending_vertices).size());

  EXPECT_EQ(starting_triangles.size(), ending_triangles.size());
  EXPECT_EQ(0, subtract(starting_triangles, ending_triangles).size());

  EXPECT_EQ(starting_sets.size(), ending_sets.size());
  EXPECT_EQ(0, subtract(starting_sets, ending_sets).size());

}
