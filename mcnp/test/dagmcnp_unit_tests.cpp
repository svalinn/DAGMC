// FluDAG/src/test/test_FlukaFuncs.cpp

#include <gtest/gtest.h>
#include "mcnp_funcs.h"

#include <cmath>
#include <cassert>

std::string test_file = "test_geom_legacy.h5m";

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
};

// Test setup outcomes
TEST_F(DAGMCNP5Test, dagmclcad_test)
{
  setup_problem();
}

// Test setup outcomes
TEST_F(DAGMCNP5Test, dagmcinit_test)
{
  std::string dagfile = test_file; 
  char* dfile = &test_file[0];
  std::string lfile = "lcad";
  char* lcadfile = &lfile[0];
  int llen = 4;
  dagmcwritemcnp_(dfile, lcadfile,&llen);
}
