// FluDAG/src/test/test_FlukaFuncs.cpp

#include <gtest/gtest.h>

#include "DagMC.hpp"
#include "moab/Interface.hpp"
#include "dagmcmetadata.hpp"

#include <cmath>
#include <cassert>

// dagmc instance
moab::DagMC *DAG = new moab::DagMC();

// metadata instance
dagmcMetaData *dgm;

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class DagmcMetadataTest : public ::testing::Test
{
 protected:

  // initalize variables for each test
  virtual void SetUp() {
    // Default h5m file for testing
    std::string infile = "test_dagmciface.h5m";

    rloadval = DAG->load_file(infile.c_str());
    assert(rloadval == moab::MB_SUCCESS);

    // DAG call to initialize geometry
    rval = DAG->init_OBBTree();
    assert (rval == moab::MB_SUCCESS);
  }

 protected:

  moab::ErrorCode rloadval;
  moab::ErrorCode rval;
};

//---------------------------------------------------------------------------//
// Test setup outcomes
TEST_F(DagmcMetadataTest, SetUp)
{
  EXPECT_EQ(moab::MB_SUCCESS, rloadval);
  // DAG call to initialize geometry
  EXPECT_EQ(moab::MB_SUCCESS, rval);
}

//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that all volumes have succesfully 
// been assigned and succesfully retreved from the metadata class
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest,TestMatAssigns)
{
  // new metadata instance
  dgm = new dagmcMetaData(DAG);
  
  // process 
  dgm->load_property_data();
  
  std::string base_property = "Hydrogen";
  std::string impl_comp_prop = "Vacuum";
  
  int num_vols = DAG->num_entities(3);
  for ( int i = 1 ; i <= num_vols ; i++ ) {
    moab::EntityHandle eh = DAG->entity_by_index(3,i);
    std::string mat_prop = dgm->get_volume_property("material",eh);

    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(mat_prop,base_property);
    else
      EXPECT_EQ(mat_prop,impl_comp_prop);

    mat_prop = dgm->get_volume_property("material",i,true);

    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(mat_prop,base_property);
    else
      EXPECT_EQ(mat_prop,impl_comp_prop);

    int cellid = DAG->id_by_index( 3, i );
    mat_prop = dgm->get_volume_property("material",i,true);

    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(mat_prop,base_property);
    else
      EXPECT_EQ(mat_prop,impl_comp_prop);
  }
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that all densities have succesfully 
// been assigned and succesfully retreved from the metadata class
// in this test there was no density data assigned, so it should be ""
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest,TestDensityAssigns)
{
  // new metadata instance
  dgm = new dagmcMetaData(DAG);
  
  // process 
  dgm->load_property_data();
  
  std::string base_property = "";
  
  int num_vols = DAG->num_entities(3);
  for ( int i = 1 ; i <= num_vols ; i++ ) {
    moab::EntityHandle eh = DAG->entity_by_index(3,i);
    std::string mat_prop = dgm->get_volume_property("density",eh);
    EXPECT_EQ(mat_prop,base_property);

    mat_prop = dgm->get_volume_property("density",i,true);
    EXPECT_EQ(mat_prop,base_property);
    
    int cellid = DAG->id_by_index( 3, i );
    mat_prop = dgm->get_volume_property("density",i,true);
    EXPECT_EQ(mat_prop,base_property);
  }
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that all volumes have succesfully 
// been assigned and succesfully retreved from the metadata class - this test
// is asserting that we have the full uwuw form for the uwuw map, i.e.
// mat:+material_name + / rho:density
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest,TestMatDensityAssigns)
{
  // new metadata instance
  dgm = new dagmcMetaData(DAG);
  
  // process 
  dgm->load_property_data();
  
  std::string base_property = "mat:Hydrogen";
  std::string impl_comp_prop = "mat:Vacuum";
  
  int num_vols = DAG->num_entities(3);
  for ( int i = 1 ; i <= num_vols ; i++ ) {
    moab::EntityHandle eh = DAG->entity_by_index(3,i);
    std::string mat_prop = dgm->get_volume_property("material_density",eh);

    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(mat_prop,base_property);
    else
      EXPECT_EQ(mat_prop,impl_comp_prop);

    mat_prop = dgm->get_volume_property("material_density",i,true);

    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(mat_prop,base_property);
    else
      EXPECT_EQ(mat_prop,impl_comp_prop);

    int cellid = DAG->id_by_index( 3, i );
    mat_prop = dgm->get_volume_property("material_density",i,true);

    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(mat_prop,base_property);
    else
      EXPECT_EQ(mat_prop,impl_comp_prop);
  }
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that all volumes have succesfully 
// been assigned and succesfully retreved from the metadata class - this test
// is asserting that we have set and correctly retrived importance data
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest,TestImportanceAssigns)
{
  // new metadata instance
  dgm = new dagmcMetaData(DAG);
  
  // process 
  dgm->load_property_data();
  
  std::string base_property = "Neutron/1.0";
  std::string impl_comp_prop = "";
  
  int num_vols = DAG->num_entities(3);
  for ( int i = 1 ; i <= num_vols ; i++ ) {
    moab::EntityHandle eh = DAG->entity_by_index(3,i);

    std::string  mat_prop = dgm->get_volume_property("importance",i,true);
    std::vector<std::string> imps = dgm->unpack_string(mat_prop,"|");
    if (!DAG->is_implicit_complement(eh))
      EXPECT_EQ(imps[0],base_property);
    else
      EXPECT_EQ(imps[0],impl_comp_prop);
  }
}
//---------------------------------------------------------------------------//
// FIXTURE-BASED TESTS: Tests to make sure that all surfaces have succesfully 
// been assigned and succesfully retreved from the dataset, specifically querying
// the boundary condition case
//---------------------------------------------------------------------------//
TEST_F(DagmcMetadataTest,TestBoundaryAssigns)
{
  // new metadata instance
  dgm = new dagmcMetaData(DAG);
  
  // process 
  dgm->load_property_data();
  
  std::string base_property = "Reflecting";
  std::string impl_comp_prop = "";
  
  int num_surfs = DAG->num_entities(2);
  int tmp[] = {1,2,3,5,6,7,8,9,11,13,14,15,16,17};
  std::vector<int> surf_ids( tmp, tmp+14 );
  for ( int i = 0 ; i < surf_ids.size(); i++ ) {
    int id = surf_ids[i];
    std::string bound_prop = dgm->get_surface_property("boundary",id,false);
    moab::EntityHandle eh = DAG->entity_by_id(2,id);
    std::string bound_prop2 = dgm->get_surface_property("boundary",eh);

    EXPECT_EQ(bound_prop,base_property);
    EXPECT_EQ(bound_prop2,base_property);

  }
}
// assert some behaviors
