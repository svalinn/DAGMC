#include <gtest/gtest.h>

#include "DagMC.hpp"
#include "moab/Interface.hpp"
#include "dagmcmetadata.hpp"

#include <cmath>
#include <cassert>

// dagmc instance
moab::DagMC* DAG;

// metadata instance
dagmcMetaData* dgm;

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class DagmcMetadataTest : public ::testing::Test {
 protected:

  // initalize variables for each test
  virtual void SetUp() {
    // Default h5m file for testing
    std::string infile = "test_dagmciface.h5m";

    DAG = new moab::DagMC();

    rloadval = DAG->load_file(infile.c_str());
    assert(rloadval == moab::MB_SUCCESS);

    // DAG call to initialize geometry
    rval = DAG->init_OBBTree();
    assert(rval == moab::MB_SUCCESS);
  }

  virtual void TearDown() {
    delete DAG;
    delete dgm;
  }

 protected:

  moab::ErrorCode rloadval;
  moab::ErrorCode rval;
};

//---------------------------------------------------------------------------//
// Test setup outcomes
TEST_F(DagmcMetadataTest, SetUp) {
  EXPECT_EQ(moab::MB_SUCCESS, rloadval);
  // DAG call to initialize geometry
  EXPECT_EQ(moab::MB_SUCCESS, rval);
}

//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that all volumes have succesfully
// been assigned and succesfully retreved from the metadata class
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest, TestMatAssigns) {
  // new metadata instance
  dgm = new dagmcMetaData(DAG);

  // process
  dgm->load_property_data();

  std::string base_property = "Hydrogen";
  std::string impl_comp_prop = "Vacuum";

  int num_vols = DAG->num_entities(3);
  for (int i = 1 ; i <= num_vols ; i++) {
    moab::EntityHandle eh = DAG->entity_by_index(3, i);
    std::string mat_prop = dgm->get_volume_property("material", eh);

    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(mat_prop, base_property);
    else
      EXPECT_EQ(mat_prop, impl_comp_prop);

    mat_prop = dgm->get_volume_property("material", i, true);

    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(mat_prop, base_property);
    else
      EXPECT_EQ(mat_prop, impl_comp_prop);

  }
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that all densities have succesfully
// been assigned and succesfully retreved from the metadata class
// in this test there was no density data assigned, so it should be ""
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest, TestDensityAssigns) {
  // new metadata instance
  dgm = new dagmcMetaData(DAG);

  // process
  dgm->load_property_data();

  std::string base_property = "";

  int num_vols = DAG->num_entities(3);
  for (int i = 1 ; i <= num_vols ; i++) {
    moab::EntityHandle eh = DAG->entity_by_index(3, i);
    std::string mat_prop = dgm->get_volume_property("density", eh);
    EXPECT_EQ(mat_prop, base_property);

    mat_prop = dgm->get_volume_property("density", i, true);
    EXPECT_EQ(mat_prop, base_property);

  }
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that all volumes have succesfully
// been assigned and succesfully retreved from the metadata class - this test
// is asserting that we have the full uwuw form for the uwuw map, i.e.
// mat:+material_name + / rho:density
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest, TestMatDensityAssigns) {
  // new metadata instance
  dgm = new dagmcMetaData(DAG);

  // process
  dgm->load_property_data();

  std::string base_property = "mat:Hydrogen";
  std::string impl_comp_prop = "mat:Vacuum";

  int num_vols = DAG->num_entities(3);
  for (int i = 1 ; i <= num_vols ; i++) {
    moab::EntityHandle eh = DAG->entity_by_index(3, i);
    std::string mat_prop = dgm->get_volume_property("material_density", eh);

    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(mat_prop, base_property);
    else
      EXPECT_EQ(mat_prop, impl_comp_prop);

    mat_prop = dgm->get_volume_property("material_density", i, true);

    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(mat_prop, base_property);
    else
      EXPECT_EQ(mat_prop, impl_comp_prop);

  }
}

TEST_F(DagmcMetadataTest, TestUnpackString) {
  // new metadata instance
  dgm = new dagmcMetaData(DAG);

  // process
  dgm->load_property_data();

  std::string neutron_property = "Neutron/1.0";
  std::string photon_property = "Photon/1.0";

  std::string  mat_prop = dgm->get_volume_property("importance", 1, true);
  std::vector<std::string> imps = dgm->unpack_string(mat_prop, "|");
  std::cout << imps[0] << std::endl;
  std::cout << imps[1] << std::endl;
  EXPECT_EQ(imps.size(), 2);
  EXPECT_EQ(imps[0], neutron_property);
  EXPECT_EQ(imps[1], photon_property);

}

//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that all volumes have succesfully
// been assigned and succesfully retreved from the metadata class - this test
// is asserting that we have set and correctly retrived importance data
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest, TestImportanceAssigns) {
  // new metadata instance
  dgm = new dagmcMetaData(DAG);

  // process
  dgm->load_property_data();

  std::string base_property = "Neutron/1.0";
  std::string impl_comp_prop = "";

  int num_vols = DAG->num_entities(3);
  for (int i = 1 ; i <= num_vols ; i++) {
    moab::EntityHandle eh = DAG->entity_by_index(3, i);

    std::string  mat_prop = dgm->get_volume_property("importance", i, true);
    std::vector<std::string> imps = dgm->unpack_string(mat_prop, "|");
    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(imps[0], base_property);
    else
      EXPECT_EQ(imps[0], impl_comp_prop);
  }
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that all surfaces have succesfully
// been assigned and succesfully retreved from the dataset, specifically querying
// the boundary condition case
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest, TestBoundaryAssigns) {
  // new metadata instance
  dgm = new dagmcMetaData(DAG);

  // process
  dgm->load_property_data();

  std::string base_property = "Reflecting";
  std::string impl_comp_prop = "";

  int num_surfs = DAG->num_entities(2);
  int tmp[] = {1, 2, 3, 5, 6, 7, 8, 9, 11, 13, 14, 15, 16, 17};
  std::vector<int> surf_ids(tmp, tmp + 14);
  for (int i = 0 ; i < surf_ids.size(); i++) {
    int id = surf_ids[i];
    std::string bound_prop = dgm->get_surface_property("boundary", id, false);
    moab::EntityHandle eh = DAG->entity_by_id(2, id);
    std::string bound_prop2 = dgm->get_surface_property("boundary", eh);

    EXPECT_EQ(bound_prop, base_property);
    EXPECT_EQ(bound_prop2, base_property);

  }
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that all surfaces have succesfully
// been assigned and succesfully retreved from the dataset, specifically querying
// the boundary condition case
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest, TestTallyAssigns) {
  // new metadata instance
  dgm = new dagmcMetaData(DAG);

  // process
  dgm->load_property_data();

  std::string base_property = "Neutron/Flux";

  int num_vols = DAG->num_entities(3);
  int tmp[] = {1, 2, 3};
  std::vector<int> vol_ids(tmp, tmp + 3);
  for (int i = 0 ; i < vol_ids.size(); i++) {
    int id = vol_ids[i];
    std::string vol_prop = dgm->get_volume_property("tally", id, false);
    std::vector<std::string> tally_props = dgm->unpack_string(vol_prop);
    moab::EntityHandle eh = DAG->entity_by_id(3, id);
    std::string vol_prop2 = dgm->get_surface_property("tally", eh);
    std::vector<std::string> tally_props2 = dgm->unpack_string(vol_prop2);

    EXPECT_EQ(tally_props[0], base_property);
    EXPECT_EQ(tally_props2[0], base_property);
  }
}

//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that the return_property function
// behaves as it is intended
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest, TestReturnProperty) {
  // new metadata instance
  dgm = new dagmcMetaData(DAG);

  std::string return_string = "";
  return_string = dgm->return_property("mat:Steel", "mat", ":", false);
  EXPECT_EQ(return_string, "mat:Steel");
  return_string = dgm->return_property("mat:Steel", "rho", ":", false);
  EXPECT_EQ(return_string, "");
  return_string = dgm->return_property("mat:Steel/rho:1.8", "mat", ":", false);
  EXPECT_EQ(return_string, "mat:Steel");
  return_string = dgm->return_property("mat:Steel/rho:1.8", "rho", ":", false);
  EXPECT_EQ(return_string, "rho:1.8");

  return_string = dgm->return_property("mat:Steel", "mat", ":", true);
  EXPECT_EQ(return_string, "Steel");
  return_string = dgm->return_property("mat:Steel", "rho", ":", true);
  EXPECT_EQ(return_string, "");
  return_string = dgm->return_property("mat:Steel/rho:1.8", "mat", ":", true);
  EXPECT_EQ(return_string, "Steel");
  return_string = dgm->return_property("mat:Steel/rho:1.8", "rho", ":", true);
  EXPECT_EQ(return_string, "1.8");

  return_string = dgm->return_property("mat:Steel", "mat");
  EXPECT_EQ(return_string, "Steel");
  return_string = dgm->return_property("mat:Steel", "rho");
  EXPECT_EQ(return_string, "");
  return_string = dgm->return_property("mat:Steel/rho:1.8", "mat");
  EXPECT_EQ(return_string, "Steel");
  return_string = dgm->return_property("mat:Steel/rho:1.8", "rho");
  EXPECT_EQ(return_string, "1.8");
}

//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that the split_string function
// behaves as it should
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest, TestSplitString) {
  // new metadata instance
  dgm = new dagmcMetaData(DAG);

  std::string to_split = "Neutron/1.0";
  std::pair<std::string, std::string> pair = dgm->split_string(to_split, "/");
  EXPECT_EQ(pair.first, "Neutron");
  EXPECT_EQ(pair.second, "1.0");

  // more complex example
  std::string more_complex = "|Neutron/1.0|Photon/2.0|";
  std::vector<std::string> imps = dgm->unpack_string(more_complex, "|");
  for (unsigned int i = 0 ; i < 2 ; i++) {
    std::string split = imps[i];
    pair = dgm->split_string(split, "/");
    if (i == 0) {
      EXPECT_EQ(pair.first, "Neutron");
      EXPECT_EQ(pair.second, "1.0");
    } else if (i == 1) {
      EXPECT_EQ(pair.first, "Photon");
      EXPECT_EQ(pair.second, "2.0");
    }
  }
}

// test to make sure the function try_to_make_int works
TEST_F(DagmcMetadataTest, TestTryToMakeInt) {
  // new metadata instance
  dgm = new dagmcMetaData(DAG);

  EXPECT_EQ(dgm->try_to_make_int("1"), true);
  EXPECT_EQ(dgm->try_to_make_int("1A"), false);
  EXPECT_EQ(dgm->try_to_make_int("M33"), false);

}
// assert some behaviors

class DagmcMetadataTestImplCompMat : public ::testing::Test {
 protected:

  // initalize variables for each test
  virtual void SetUp() {
    // Default h5m file for testing
    std::string infile = "test_dagmciface_impl.h5m";

    DAG = new moab::DagMC();

    rloadval = DAG->load_file(infile.c_str());
    assert(rloadval == moab::MB_SUCCESS);

    // DAG call to initialize geometry
    rval = DAG->init_OBBTree();
    assert(rval == moab::MB_SUCCESS);
  }

  virtual void TearDown() {
    //    delete dgm;
    delete DAG;
  }

 protected:

  moab::ErrorCode rloadval;
  moab::ErrorCode rval;
};

//---------------------------------------------------------------------------//
// Test setup outcomes
TEST_F(DagmcMetadataTestImplCompMat, SetUp) {
  EXPECT_EQ(moab::MB_SUCCESS, rloadval);
  // DAG call to initialize geometry
  EXPECT_EQ(moab::MB_SUCCESS, rval);
}

// make sure the the implicit complement material
// is set
TEST_F(DagmcMetadataTestImplCompMat, ImplCompMat) {
  // new metadata instance
  dgm = new dagmcMetaData(DAG);
  // process
  dgm->load_property_data();
  // loop over the volumes
  int num_vols = DAG->num_entities(3);

  std::string mat1 = "Steel";
  std::string mat_grave = "Graveyard";
  std::string mat_impl = "Steel";

  std::string mat_prop = dgm->get_volume_property("material", 1, true);
  EXPECT_EQ(mat1, mat_prop);

  std::string mat_prop2 = dgm->get_volume_property("material", 2, true);
  EXPECT_EQ(mat_grave, mat_prop2);

  std::string mat_prop3 = dgm->get_volume_property("material", 3, true);
  EXPECT_EQ(mat_impl, mat_prop3);

}
