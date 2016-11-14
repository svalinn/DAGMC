
#include "gtest/gtest.h"
#include "test_classes.hpp"

class SphereNBoxMakeWatertightTest : public MakeWatertightTest
{

 protected:
  virtual void setFilename() {
    filename = "sphere_n_box.h5m";
  };

};


TEST_F(SphereNBoxMakeWatertightTest, SphereNBoxTest )
{
  //make sure that the expected number of surfaces exist
  int dim = 2, expected_num_surfs = 7;
  moab::ErrorCode rval;
  rval = check_num_ents(dim, expected_num_surfs);
  EXPECT_TRUE(moab::MB_SUCCESS == rval);
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
  rval = check_num_ents(dim, expected_num_surfs);
  EXPECT_TRUE(moab::MB_SUCCESS == rval);
}
