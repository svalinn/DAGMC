#include <gtest/gtest.h>

#include <filesystem>
#include <iostream>
#include <memory>

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

// This test exists to ensure that the IPC is correctly positioned in the model
// See GitHub issue #934 for more information
// The test loads a known, working DAGMC file, adds a bunch of meshsets to it,
// and then writes a temporary file. The temporary file is then read back in
// and the IPC is checked to ensure it is in the correct position. The condision
// being tested is only triggered during a file read when EntityHandle sequences
// are generated and accessed internally and cannot be reproduced
// using the external-facing MOAB API.
TEST_F(DagmcIPCPositionTest, dagmc_implicit_complement_position_test) {
  // create a DAGMC instance
  std::unique_ptr<DagMC> dagmc = std::make_unique<DagMC>();

  // load the geometry test file
  ErrorCode rval = dagmc->load_file(simple_file.c_str());

  // add one million mesh sets to this file to ensure the threshold of
  // meshsets in a file is surpassed to trigger this condition
  // NOTE: far fewer meshsets are needed to trigger this condition, but
  //       one million was chosen to be conservative should the internal
  //       allocation size in MOAB change again in the future
  EntityHandle tmp;
  for (int i = 0; i < 1E6; i++) {
    rval = dagmc->moab_instance()->create_meshset(MESHSET_SET, tmp);
    ASSERT_EQ(rval, MB_SUCCESS);
  }

  // write a temprary file
  dagmc->moab_instance()->write_file("tmp.h5m");

  // create a second DAGMC instance and read the temporary file
  std::unique_ptr<DagMC> dagmc2 = std::make_unique<DagMC>();
  dagmc2->load_file("tmp.h5m");

  // perform necessary DAGMC setup
  dagmc2->geom_tool()->find_geomsets();
  dagmc2->setup_impl_compl();
  dagmc2->setup_indices();

  // get the implicit complement handle, it should exist after the explicit call
  // for its creation above
  EntityHandle ipc = 0;
  rval = dagmc2->geom_tool()->get_implicit_complement(ipc);
  ASSERT_EQ(rval, MB_SUCCESS);
  ASSERT_NE(ipc, 0);

  // there are 3 volumes in the original model and we've added the implicit
  // complement
  int num_vols = dagmc2->num_entities(3);
  ASSERT_EQ(num_vols, 4);

  // Make sure the IPC handle index is highest
  // Reminder: DAGMC indieces are 1-based
  std::cout << "IPC index: " << dagmc2->index_by_handle(ipc) << std::endl;
  ASSERT_EQ(dagmc2->index_by_handle(ipc), 4);
}