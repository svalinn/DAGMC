
#include <memory>

#include "gtest/gtest.h"
#include "moab/Core.hpp"

using namespace moab;

class OverlapTest : public ::testing::Test {
 protected:
  virtual void SetUp() override;
  virtual void TearDown() override;
  virtual void SetFilename(){};

  std::string filename;

 public:
  std::shared_ptr<Interface> MBI;
};