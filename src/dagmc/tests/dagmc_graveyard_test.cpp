#include <gtest/gtest.h>

#include <iostream>
#include <memory>

#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"

using namespace moab;

using moab::DagMC;

static std::string simple_file = "test_dagmc.h5m";
static std::string trelis_file = "pincell.h5m";

class DagmcGraveyardTest : public ::testing::Test {
 protected:
  virtual void SetUp() {}
  virtual void TearDown() {}
};

TEST_F(DagmcGraveyardTest, dagmc_graveyard_simple_test) {
  // create a DAGMC instance
  std::unique_ptr<DagMC> DAG(new DagMC());
  // load the test geometry file
  ErrorCode rval = DAG->load_file(simple_file.c_str());  // open the Dag file
  EXPECT_EQ(MB_SUCCESS, rval);

  rval = DAG->init_OBBTree();
  EXPECT_EQ(MB_SUCCESS, rval);

  // collect starting sets to make sure we end up with the same thing at the end
  Range starting_vertices;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBVERTEX,
                                                    starting_vertices);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range starting_triangles;
  rval =
      DAG->moab_instance()->get_entities_by_type(0, MBTRI, starting_triangles);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range starting_sets;
  rval =
      DAG->moab_instance()->get_entities_by_type(0, MBENTITYSET, starting_sets);
  EXPECT_EQ(MB_SUCCESS, rval);

  // there is no graveyard present to start in this model
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
  rval =
      DAG->moab_instance()->get_entities_by_type(0, MBVERTEX, ending_vertices);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range ending_triangles;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBTRI, ending_triangles);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range ending_sets;
  rval =
      DAG->moab_instance()->get_entities_by_type(0, MBENTITYSET, ending_sets);
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_EQ(starting_vertices.size(), ending_vertices.size());
  EXPECT_EQ(0, subtract(starting_vertices, ending_vertices).size());

  EXPECT_EQ(starting_triangles.size(), ending_triangles.size());
  EXPECT_EQ(0, subtract(starting_triangles, ending_triangles).size());

  EXPECT_EQ(starting_sets.size(), ending_sets.size());
  EXPECT_EQ(0, subtract(starting_sets, ending_sets).size());
}

TEST_F(DagmcGraveyardTest, dagmc_graveyard_test_trelis_file) {
  // create DAGMC instance
  std::unique_ptr<DagMC> DAG(new DagMC());
  // load the trelis test file
  ErrorCode rval = DAG->load_file(trelis_file.c_str());
  EXPECT_EQ(MB_SUCCESS, rval);

  rval = DAG->init_OBBTree();
  EXPECT_EQ(MB_SUCCESS, rval);

  // collect starting sets to make sure we end up with the same thing at the end
  Range starting_vertices;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBVERTEX,
                                                    starting_vertices);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range starting_triangles;
  rval =
      DAG->moab_instance()->get_entities_by_type(0, MBTRI, starting_triangles);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range starting_sets;
  rval =
      DAG->moab_instance()->get_entities_by_type(0, MBENTITYSET, starting_sets);
  EXPECT_EQ(MB_SUCCESS, rval);

  // a graveyard is already present in this model
  EXPECT_TRUE(DAG->has_graveyard());

  int n_groups = DAG->num_entities(4);
  int n_vols = DAG->num_entities(3);
  int n_surfs = DAG->num_entities(2);
  // the DagMC class only tracks relevant geometry sets (surfaces and volumes)
  int n_curves = DAG->geom_tool()->num_ents_of_dim(1);
  int n_geom_verts = DAG->geom_tool()->num_ents_of_dim(0);

  // remove original graveyard
  rval = DAG->remove_graveyard();
  EXPECT_EQ(MB_SUCCESS, rval);

  // geometric sets removed:
  // (geometric vertices - 16, curves - 24, surfaces - 12, volumes - 1, groups -
  // 1)
  EXPECT_EQ(16, n_geom_verts - DAG->geom_tool()->num_ents_of_dim(0));
  EXPECT_EQ(24, n_curves - DAG->geom_tool()->num_ents_of_dim(1));
  EXPECT_EQ(12, n_surfs - DAG->num_entities(2));
  EXPECT_EQ(1, n_vols - DAG->num_entities(3));
  EXPECT_EQ(1, n_groups - DAG->num_entities(4));

  DAG->moab_instance()->write_file("test.h5m");
  // set of vertices, triangles, and sets without the original graveyard volume
  Range model_vertices;
  rval =
      DAG->moab_instance()->get_entities_by_type(0, MBVERTEX, model_vertices);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range model_triangles;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBTRI, model_triangles);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range model_sets;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBENTITYSET, model_sets);
  EXPECT_EQ(MB_SUCCESS, rval);

  // make sure the right number of mesh elements were removed
  // sixteen vertices removed (corners of two cuboid volumes)
  EXPECT_EQ(starting_vertices.size() - 16, model_vertices.size());
  // twenty-four triangles removed (two per face for two cuboid volumes)
  EXPECT_EQ(starting_triangles.size() - 24, model_triangles.size());

  // update number of surfaces and volumes before creating the graveyard
  n_vols = DAG->num_entities(3);
  n_surfs = DAG->num_entities(2);

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
  rval =
      DAG->moab_instance()->get_entities_by_type(0, MBVERTEX, ending_vertices);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range ending_triangles;
  rval = DAG->moab_instance()->get_entities_by_type(0, MBTRI, ending_triangles);
  EXPECT_EQ(MB_SUCCESS, rval);

  Range ending_sets;
  rval =
      DAG->moab_instance()->get_entities_by_type(0, MBENTITYSET, ending_sets);
  EXPECT_EQ(MB_SUCCESS, rval);

  EXPECT_EQ(model_vertices.size(), ending_vertices.size());
  EXPECT_EQ(0, subtract(model_vertices, ending_vertices).size());

  EXPECT_EQ(model_triangles.size(), ending_triangles.size());
  EXPECT_EQ(0, subtract(model_triangles, ending_triangles).size());

  EXPECT_EQ(model_sets.size(), ending_sets.size());
  EXPECT_EQ(0, subtract(model_sets, ending_sets).size());
}
