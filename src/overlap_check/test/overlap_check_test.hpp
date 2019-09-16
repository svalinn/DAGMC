
#include "gtest/gtest.h"

using namespace moab;

class OverlapTest : public::testing::Test {
protected:
  virtual void SetUp();
  virtual void TearDown();
  virtual void SetFilename() {};

  std::string filename;

public:
  std::shared_ptr<Interface> MBI;
};
