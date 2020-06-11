
#include "gtest/gtest.h"

using namespace moab;

class OverlapTest : public::testing::Test {
 protected:
  virtual void SetUp() override;
  virtual void TearDown() override;
  virtual void SetFilename() override {};

  std::string filename;

 public:
  std::shared_ptr<Interface> MBI;
};
