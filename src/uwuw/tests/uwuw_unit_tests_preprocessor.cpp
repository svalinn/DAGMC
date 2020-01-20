
#include <gtest/gtest.h>
#include <cmath>
#include <cassert>

#include <iostream>
#include <unistd.h>
#include <stdio.h>

#include "uwuw_preprocessor.hpp"
#include "pyne.h"

namespace {

class UWUWTest : public ::testing::Test {
 protected:

  UWUWTest() {}
  virtual ~UWUWTest() {}

  virtual void SetUp() {
  }

};


/*
 * Empty common setup function
 */
TEST_F(UWUWTest, SetUp) {
}

/*
 * Test to make sure the name is only uppercased
 */
TEST_F(UWUWTest, name8) {
  // make new class
  name_concatenator* ncr = new name_concatenator();

  std::string name = "SpecialSteel";
  EXPECT_EQ(ncr->make_name_8bytes(name), "SPECIALS");
  return;
}

// test to make sure the space is removed
TEST_F(UWUWTest, name8RemoveSpace) {
  // make new class
  name_concatenator* ncr = new name_concatenator();

  std::string name = "Special Steel";
  EXPECT_EQ(ncr->make_name_8bytes(name), "SPECIALS");
  return;
}

// test to make sure the , is removed
TEST_F(UWUWTest, name8RemoveCommaAndSpace) {
  // make new class
  name_concatenator* ncr = new name_concatenator();

  std::string name = "Special, Steel";
  EXPECT_EQ(ncr->make_name_8bytes(name), "SPECIALS");
  return;
}

// test to make sure the number is incremented
TEST_F(UWUWTest, name8IncrementID) {
  // make new class
  name_concatenator* ncr = new name_concatenator();

  std::string name = "Special, Steel";
  EXPECT_EQ(ncr->make_name_8bytes(name), "SPECIALS");
  EXPECT_EQ(ncr->make_name_8bytes(name), "SPECIAL1");
  return;
}

// test to make sure the number is incremented
TEST_F(UWUWTest, name8lessthan8) {
  // make new class
  name_concatenator* ncr = new name_concatenator();

  std::string name = "Air, Dry";
  EXPECT_EQ(ncr->make_name_8bytes(name), "  AIRDRY");
  return;
}

// test to make sure the number is incremented but that
// we preserve whitespace
TEST_F(UWUWTest, name8lessthan8preserve) {
  // make new class
  name_concatenator* ncr = new name_concatenator();

  std::string name = "Air, Dry";
  EXPECT_EQ(ncr->make_name_8bytes(name), "  AIRDRY");
  EXPECT_EQ(ncr->make_name_8bytes(name), " AIRDRY1");

  return;
}

// test to make sure the number is incremented but that
// we preserve whitespace and the number is good
TEST_F(UWUWTest, name8lessthan8preserve100) {
  // make new class
  name_concatenator* ncr = new name_concatenator();

  std::string name = "Air, Dry";
  EXPECT_EQ(ncr->make_name_8bytes(name), "  AIRDRY");

  for (int i = 1 ; i < 10 ; i++) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name), " AIRDRY" + str);
  }
  for (int i = 10 ; i < 100 ; i++) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name), "AIRDRY" + str);
  }
  for (int i = 100 ; i < 999 ; i++) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name), "AIRDR" + str);
  }

  return;
}

// test to make sure the number is incremented but that
// we preserve whitespace and the number is good
TEST_F(UWUWTest, name8lessthan8preserve1000Steel) {
  // make new class
  name_concatenator* ncr = new name_concatenator();

  std::string name = "Special Steel";
  EXPECT_EQ(ncr->make_name_8bytes(name), "SPECIALS");

  for (int i = 1 ; i < 10 ; i++) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name), "SPECIAL" + str);
  }
  for (int i = 10 ; i < 100 ; i++) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name), "SPECIA" + str);
  }
  for (int i = 100 ; i < 999 ; i++) {
    std::stringstream ss;
    ss << i;

    std::string str;
    ss >> str;
    EXPECT_EQ(ncr->make_name_8bytes(name), "SPECI" + str);
  }

  return;
}

// test to ensure that when a material is taken
// from the library file and inserted into uwuw material
// we bring all the metadat with us
TEST_F(UWUWTest, materialMetadata) {
  std::string lib_file = "mat_lib.h5";
  std::string dag_file = "dag_file.h5m";
  std::string out_file = "intermediate.h5";
  std::string matlib_hdf5_path = "/materials";
  bool verbose = false;
  bool fatal_errors = false;
  // make new preprocessor
  uwuw_preprocessor* uwuw_preproc = new uwuw_preprocessor(lib_file, dag_file,
                                                          out_file, matlib_hdf5_path,
                                                          verbose, fatal_errors);
  // load the geometry
  // process materials
  uwuw_preproc->process_materials();
  // write the new material library
  uwuw_preproc->write_uwuw_materials();

  // now read in the material library
  UWUW* uwuw = new UWUW(out_file);
  std::map<std::string, pyne::Material> mat_lib = uwuw->material_library;

  // pull out the only material
  pyne::Material mat = mat_lib["mat:CentreStack"];

  EXPECT_EQ(mat.metadata["name"].asString(), "mat:CentreStack");
  EXPECT_EQ(mat.metadata["fluka_name"].asString(), "CENTREST");
  EXPECT_EQ(mat.metadata["mat_number"].asInt(), 1);
  EXPECT_EQ(mat.metadata["special_tag"].asString(), "this is a test tag");

  delete uwuw;
  delete uwuw_preproc;
  std::remove(out_file.c_str());
}
};


namespace {

class UWUWPreprocTest : public ::testing::Test {
 protected:
  uwuw_preprocessor* uwuw_preproc;

  UWUWPreprocTest() {}
  virtual ~UWUWPreprocTest() {}

  virtual void SetUp() {
    // make new preprocessor
    uwuw_preproc = new uwuw_preprocessor("mat_lib.h5",
                                         "dag_file.h5m",
                                         "output_test_file.h5",
                                         "/materials",
                                         true,
                                         false);

    // process the materials
    uwuw_preproc->process_materials();

    // process the tallies
    uwuw_preproc->process_tallies();

  }

  virtual void TearDown() {
    delete uwuw_preproc;
  }
};

/*
 * Test to make sure the name is only uppercased
 */
TEST_F(UWUWPreprocTest, testMaterials) {
  // write the material data
  uwuw_preproc->write_uwuw_materials();

  // write the tally data
  uwuw_preproc->write_uwuw_tallies();

  std::remove("output_test_file.h5m");
}

/*
 * Test to make sure the name is only uppercased
 */
TEST_F(UWUWPreprocTest, testTallies) {
  // write the material data
  uwuw_preproc->write_uwuw_materials();

  // write the tally data
  uwuw_preproc->write_uwuw_tallies();

  std::remove("output_test_file.h5m");
}
};
