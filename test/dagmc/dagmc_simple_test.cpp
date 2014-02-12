#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"

#include "DagMC.hpp"

using namespace moab;

using moab::DagMC;

#define DAG DagMC::instance()

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)

#ifdef MESHDIR
static const char input_file[] = STRINGIFY(MESHDIR) "/dagmc/test_geom.h5m";
#else
static const char input_file[] = STRINGIFY(MESHDIR) "/dagmc/test_geom.h5m";
#endif


void dagmc_load_file() 
{
  ErrorCode rval = DAG->load_file(input_file); // open the Dag file
  CHECK_ERR(rval);
}

void dagmc_build_obb() 
{
  ErrorCode rval = DAG->init_OBBTree();
  CHECK_ERR(rval);
}

void dagmc_num_vols()
{
  int num_vols = DAG->num_entities(3); 
  CHECK_EQUAL(2,num_vols);
}

void dagmc_entity_handle()
{
  int num_vols = DAG->num_entities(3); 
  EntityHandle vol;
  EntityHandle volume = 12682136550675316765;
  for ( int i = 0 ; i < num_vols ; i++ )
    {
      vol = DAG->entity_by_index(3,i);
    }
  CHECK_EQUAL(volume,vol);
}

void dagmc_point_in()
{
  int result = 0;
  double xyz[3]={0.0,0.0,0.0};
  EntityHandle volume = 12682136550675316765;
  ErrorCode rval = DAG->point_in_volume(volume,xyz,result);
  CHECK_EQUAL(1,result);
}

void dagmc_rayfire()
{
  const double eps = 1e-6; // epsilon for test, faceting tol?

  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3]={0.0,0.0,0.0};
  double dir[3]={0.0,0.0,1.0}; 
  EntityHandle next_surf;
  double next_surf_dist;
  EntityHandle volume = 12682136550675316765;

  ErrorCode rval = DAG->ray_fire(volume,xyz,dir,next_surf,next_surf_dist);
  CHECK_REAL_EQUAL(5.0,next_surf_dist,eps); 
}

void dagmc_closest_to()
{
  const double eps = 1e-6; // epsilon for test, faceting tol?

  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3]={-6.0,0.0,0.0};
  double distance; // distance from point to nearest surface
  EntityHandle volume = 12682136550675316765;

  ErrorCode rval = DAG->closest_to_location(volume,xyz,distance);
  // distance should be 1.0 cm
  CHECK_REAL_EQUAL(1.0,distance,eps); 
}

void dagmc_test_boundary()
{
  EntityHandle volume = 12682136550675316765;
  EntityHandle surf   = 12682136550675316759;
  double xyz[3]={0.0,0.0,5.0};
  double dir[3]={0.0,0.0,1.0}; 
  int result;
  
  ErrorCode rval = DAG->test_volume_boundary(volume,surf,xyz,dir,result);
  // check ray leaving volume
  CHECK_EQUAL(0,result);  
}
  
int main(int /* argc */, char** /* argv */)
{
  int result = 0;
  result += RUN_TEST( dagmc_load_file );     // test ray fire
  result += RUN_TEST( dagmc_build_obb );     // build the obb
  result += RUN_TEST( dagmc_num_vols  );     // make sure the num of vols correct
  result += RUN_TEST( dagmc_entity_handle);  // check the entity handle correct
  result += RUN_TEST( dagmc_point_in);       // check entity by point
  result += RUN_TEST( dagmc_rayfire ) ;      // ensure ray fire distance is correct
  result += RUN_TEST( dagmc_closest_to );    // check the distance to surface nearest point
  result += RUN_TEST( dagmc_test_boundary ); // check particle entering leaving

  return result;
}
