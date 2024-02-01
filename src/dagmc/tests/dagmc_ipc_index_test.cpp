#include <gtest/gtest.h>

#include <iostream>
#include <memory>
#include <filesystem>

#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"

static std::string simple_file = "test_dagmc.h5m";

class DagmcIPCPositionTest : public ::testing::Test {
 protected:
  virtual void SetUp() {}
  virtual void TearDown() {
    if (std::filesystem::exists("tmp.h5m")) {
      std::filesystem::remove("tmp.h5m");
    }
  }
};

using namespace moab;


TEST_F(DagmcIPCPositionTest, dagmc_implicit_complement_position_test) {
  // create a DAGMC instance
  std::unique_ptr<DagMC> dagmc = std::make_unique<DagMC>();

  // load the geometry test file
  ErrorCode rval = dagmc->load_file(simple_file.c_str());

  // add a bunch of meshsets to this file
  EntityHandle tmp;
  for (int i = 0; i < 1E6; i++) {
    rval = dagmc->moab_instance()->create_meshset(MESHSET_SET, tmp);
    ASSERT_EQ(rval, MB_SUCCESS);
  }

  dagmc->moab_instance()->write_file("tmp.h5m");

  // create a second DAGMC instance nd read the new file
  std::unique_ptr<DagMC> dagmc2 = std::make_unique<DagMC>();

  dagmc2->load_file("tmp.h5m");
  dagmc2->geom_tool()->find_geomsets();
  dagmc2->setup_impl_compl();
  dagmc2->setup_indices();

  EntityHandle ipc = 0;
  rval = dagmc2->geom_tool()->get_implicit_complement(ipc);
  ASSERT_EQ(rval, MB_SUCCESS);
  ASSERT_NE(ipc, 0);

  int num_vols = dagmc2->num_entities(3);
  // 3 volumes in the original model, plus the IPC
  ASSERT_EQ(num_vols, 4);

  // Reminder: DAGMC indieces are 1-based
  std::cout << "IPC index: " << dagmc2->index_by_handle(ipc) << std::endl;
  ASSERT_EQ(dagmc2->index_by_handle(ipc), 4);
}