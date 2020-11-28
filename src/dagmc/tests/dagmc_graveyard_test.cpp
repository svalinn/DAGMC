#include <gtest/gtest.h>

#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "DagMC.hpp"

#include <iostream>
#include <memory>

using namespace moab;

using moab::DagMC;

static std::string simple_file = "test_dagmc.h5m";
static std::string trelis_file = "pincell.h5m";

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

TEST_F(DagmcGraveyardTest, dagmc_graveyard_simple_test) {
  std::unique_ptr<DagMC> DAG(new DagMC());
  ErrorCode rval = DAG->load_file(simple_file.c_str()); // open the Dag file
  EXPECT_EQ(MB_SUCCESS, rval);

  rval = DAG->init_OBBTree();
  EXPECT_EQ(MB_SUCCESS, rval);

  // collect starting sets to make sure we end up with the same thing at the end
  Range starting_vertices;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBVERTEX, starting_vertices);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range starting_triangles;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBTRI, starting_triangles);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range starting_sets;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBENTITYSET, starting_sets);
  EXPECT_EQ(MB_SUCCESS, rval);

  // there is no graveyard present in this model
  EXPECT_FALSE(DAG->has_graveyard());

  int n_vols = DAG->num_entities(3);
  int n_surfs = DAG->num_entities(2);

  rval = DAG->create_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_EQ(n_vols + 1, DAG->num_entities(3));
  EXPECT_EQ(n_surfs + 2, DAG->num_entities(2));

  EXPECT_TRUE(DAG->has_graveyard());

  rval = DAG->remove_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_EQ(n_vols, DAG->num_entities(3));
  EXPECT_EQ(n_surfs, DAG->num_entities(2));

  EXPECT_FALSE(DAG->has_graveyard());

  rval = DAG->create_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_TRUE(DAG->has_graveyard());

  // test overwrite capability
  bool overwrite_graveyard = true;
  rval = DAG->create_graveyard(true);
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_EQ(n_vols + 1, DAG->num_entities(3));
  EXPECT_EQ(n_surfs + 2, DAG->num_entities(2));

  // this should fail, graveyard already exists
  // and overwrite isn't specified
  rval = DAG->create_graveyard();
  EXPECT_EQ(MB_FAILURE, rval);

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

TEST_F(DagmcGraveyardTest, dagmc_graveyard_test_trelis_file) {
  std::unique_ptr<DagMC> DAG(new DagMC());
  ErrorCode rval = DAG->load_file(trelis_file.c_str()); // open the Dag file
  EXPECT_EQ(MB_SUCCESS, rval);

  rval = DAG->init_OBBTree();
  EXPECT_EQ(MB_SUCCESS, rval);

   // there should already be a graveyard present in this model
  EXPECT_TRUE(DAG->has_graveyard());

  int n_vols = DAG->num_entities(3);
  int n_surfs = DAG->num_entities(2);

  rval = DAG->remove_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_EQ(n_vols - 1, DAG->num_entities(3));
  // each cube face is a surface here, there should be 12 fewer
  EXPECT_EQ(n_surfs - 12, DAG->num_entities(2));

  // update number of surfaces and volumes before creating the graveyard
  n_vols = DAG->num_entities(3);
  n_surfs = DAG->num_entities(2);

  // collect starting sets to make sure we end up with the same thing at the end
  Range starting_vertices;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBVERTEX, starting_vertices);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range starting_triangles;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBTRI, starting_triangles);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range starting_sets;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBENTITYSET, starting_sets);
  EXPECT_EQ(MB_SUCCESS, rval);

  rval = DAG->create_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_EQ(n_vols + 1, DAG->num_entities(3));
  EXPECT_EQ(n_surfs + 2, DAG->num_entities(2));

  EXPECT_TRUE(DAG->has_graveyard());

  rval = DAG->remove_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_EQ(n_vols, DAG->num_entities(3));
  EXPECT_EQ(n_surfs, DAG->num_entities(2));

  EXPECT_FALSE(DAG->has_graveyard());

  rval = DAG->create_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_TRUE(DAG->has_graveyard());

  // test overwrite capability
  bool overwrite_graveyard = true;
  rval = DAG->create_graveyard(true);
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_EQ(n_vols + 1, DAG->num_entities(3));
  EXPECT_EQ(n_surfs + 2, DAG->num_entities(2));

  // this should fail, graveyard already exists
  // and overwrite isn't specified
  rval = DAG->create_graveyard();
  EXPECT_EQ(MB_FAILURE, rval);

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
