

#include "gtest/gtest.h"
#include "test_classes.hpp"

class ITERMakeWatertightRegressionTest : public MakeWatertightTest
{

 protected:
  virtual void setFilename() {
    filename = "iter_imprinted.h5m";
  };

};

class FNSFMakeWatertightRegressionTest : public MakeWatertightTest
{

 protected:
  virtual void setFilename() {
    filename = "model1_360_1.h5m";
  };

};

class BLiteMakeWatertightRegressionTest : public MakeWatertightTest
{

 protected:
  virtual void setFilename() {
    filename = "bllite30matls.h5m";
  };

};

TEST_F(ITERMakeWatertightRegressionTest, ITERRegressionTest )
{
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}

TEST_F(FNSFMakeWatertightRegressionTest, FNSFRegressionTest )
{
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}

TEST_F(BLiteMakeWatertightRegressionTest, BLiteRegressionTest )
{
  EXPECT_TRUE(seal_and_check(input_fileset, facet_tol));
}
