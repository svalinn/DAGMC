#include <gtest/gtest.h>
#include "mcnp_funcs.h"

#include <cmath>
#include <cassert>

std::string test_file = "test_geom_legacy.h5m";
std::string test_file_comp = "test_geom_legacy_comp.h5m";

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class DAGMCNP5Test : public ::testing::Test
{
 public:
  void setup_problem() {
    std::string filename = test_file;
    char* file = &filename[0];
    int len = filename.length();
    std::string facet_tol = "1.0e-4";
    char* ftol = &facet_tol[0];
    int ftol_len = facet_tol.length();
    int parallel_mode = 0;
    double dagmc_version;
    int moab_version;
    int max_pbl;

    // intialise dagmc
    dagmcinit_(file,&len,ftol,&ftol_len,&parallel_mode,
               &dagmc_version,&moab_version,&max_pbl);

  }
  void setup_problem_comp() {
    std::string filename = test_file_comp;
    char* file = &filename[0];
    int len = filename.length();
    std::string facet_tol = "1.0e-4";
    char* ftol = &facet_tol[0];
    int ftol_len = facet_tol.length();
    int parallel_mode = 0;
    double dagmc_version;
    int moab_version;
    int max_pbl;

    // intialise dagmc
    dagmcinit_(file,&len,ftol,&ftol_len,&parallel_mode,
               &dagmc_version,&moab_version,&max_pbl);

  }

};

// Test setup outcomes
TEST_F(DAGMCNP5Test, dagmclcad_test)
{
  setup_problem();
}

// Test setup outcomes
TEST_F(DAGMCNP5Test, dagmcinit_test)
{
// expected values from the lcad file // only the cells
  const char* expected[] = {"1 9 -12.0 imp:n=1 imp:p=1 ",
                            "2 9 -12.0 imp:n=1 imp:p=1 ",
                            "3 9 -12.0 imp:n=1 imp:p=1 ",
                            "4 9 -12.0 imp:n=1 imp:p=1 ",
                            "5 3 -2.0 imp:n=1 imp:p=1 ",
                            "6 3 -2.0 imp:n=1 imp:p=1 ",
                            "7 3 -2.0 imp:n=1 imp:p=1 ",
                            "8 3 -2.0 imp:n=1 imp:p=1 ",
                            "9 1 0.022 imp:n=1 imp:p=1 ",
                            "12 0  imp:n=0 imp:p=0   $ graveyard",
                            "13 0  imp:n=1 imp:p=1   $ implicit complement"
                           };
  std::vector<std::string> expected_lcad(expected,expected+11);

  std::string dagfile = test_file;
  char* dfile = &test_file[0];
  std::string lfile = "lcad";
  char* lcadfile = &lfile[0];
  int llen = 4;
  dagmcwritemcnp_(dfile, lcadfile,&llen);

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

  // clearout the dagmc instance
  dagmc_teardown_();
}

// Test setup outcomes
TEST_F(DAGMCNP5Test, dagmclcad_comp_test)
{
  setup_problem_comp();
}

// Test setup outcomes
TEST_F(DAGMCNP5Test, dagmcinit_comp_test)
{
// expected values from the lcad file // only the cells
  const char* expected[] = {"1 9 -12.0 imp:n=1 imp:p=1 ",
                            "2 9 -12.0 imp:n=1 imp:p=1 ",
                            "3 9 -12.0 imp:n=1 imp:p=1 ",
                            "4 9 -12.0 imp:n=1 imp:p=1 ",
                            "5 3 -2.0 imp:n=1 imp:p=1 ",
                            "6 3 -2.0 imp:n=1 imp:p=1 ",
                            "7 3 -2.0 imp:n=1 imp:p=1 ",
                            "8 3 -2.0 imp:n=1 imp:p=1 ",
                            "9 1 0.022 imp:n=1 imp:p=1 ",
                            "12 0  imp:n=0 imp:p=0   $ graveyard",
                            "13 2 -3.1 imp:n=1 imp:p=1   $ implicit complement"
                           };
  std::vector<std::string> expected_lcad(expected,expected+11);

  std::string dagfile = test_file;
  char* dfile = &test_file[0];
  std::string lfile = "lcad";
  char* lcadfile = &lfile[0];
  int llen = 4;
  dagmcwritemcnp_(dfile, lcadfile,&llen);

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

  // clearout the dagmc instance
  dagmc_teardown_();
}
