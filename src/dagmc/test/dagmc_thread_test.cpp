#include <gtest/gtest.h>

#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "DagMC.hpp"
#include "thread_manager.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <random>
#include <iostream>

using namespace moab;
using moab::DagMC;

DagThreadManager* DTM;

static const char input_file[] = "test_geom.h5m";

extern int num_thread_to_run;

#define TWO_PI 6.28318530718 // TODO need more sf

int num_threads = 0;

class DagmcThreadTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    num_threads = num_thread_to_run;

    // set the number of threads to be the number of procs
    DTM = new DagThreadManager(num_threads);

    // get a dagmc instance to continue setup
    moab::DagMC* DAG = DTM->get_dagmc_instance(0);
    moab::ErrorCode rloadval = DAG->load_file(input_file);
    assert(rloadval == moab::MB_SUCCESS);
    // Create the OBB
    rval = DAG->init_OBBTree();
    assert(rval == moab::MB_SUCCESS);
    //    DTM->setup_child_threads();
    DTM->initialise_child_threads();
  }
  // clearout
  //  virtual void TearDown() {
  //  delete DTM;
  //}
 protected:
  moab::ErrorCode rloadval;
  moab::ErrorCode rval;
};

// function that returns the EH of the volume the point is in
moab::EntityHandle where_am_i(moab::DagMC* DAG,
                              moab::DagMC::RayHistory* RayHist,
                              double position[3], double dir[3]) {
  int i = 0;
  int num_vols = DAG->num_entities(3);
  int inside = 0;
  for (i = 1 ; i <= num_vols ; i++) {
    moab::EntityHandle eh = DAG->entity_by_index(3, i);
    moab::ErrorCode ec = DAG->point_in_volume(eh, position, inside,
                                              dir, RayHist);
    if (inside)
      return eh;
  }
  return 0;
}

// given the surface and the current volume, where do we go next?
moab::EntityHandle next_volume(moab::DagMC* DAG, moab::EntityHandle surf,
                               moab::EntityHandle eh) {
  moab::EntityHandle next_vol = 0;
  moab::ErrorCode rval = DAG->next_vol(surf, eh, next_vol);
  assert(rval == moab::MB_SUCCESS);
  return next_vol;
}

std::vector<double> isotropic_dir(int nps) {
  // see the rn

  int stride = 123456;
  // seed the generator for each nps - unique seed
  std::mt19937 gen(nps * stride);
  std::uniform_real_distribution<> rand(0, 1);
  double r1 = rand(gen);
  double r2 = rand(gen);

  std::vector<double> dir;
  double u = 2.*r1 - 1.;
  double v = std::sqrt(1. - (u * u)) * cos(TWO_PI * r2);
  double w = std::sqrt(1. - (u * u)) * sin(TWO_PI * r2);
  dir.push_back(u);
  dir.push_back(v);
  dir.push_back(w);
  return dir;
}

// fire a ray frrom start along dir, which surf do we hit and what distance
moab::EntityHandle fire_ray(moab::DagMC* DAG, moab::DagMC::RayHistory* RayHist,
                            moab::EntityHandle eh, double start[3], double dir[3],
                            double& distance) {
  moab::EntityHandle next_surf = 0;
  moab::ErrorCode rval = DAG->ray_fire(eh, start, dir, next_surf,
                                       distance, RayHist);
  assert(rval == moab::MB_SUCCESS);
  return next_surf;
}

// mimics the MCNP transport cycle
moab::EntityHandle transport_cycle(moab::DagMC* DAG,
                     moab::DagMC::RayHistory* RayHistory,
                     int nps,
                     double start[3]) {

  std::vector<double> rand_dir = isotropic_dir(nps);
  double dir[3]; dir[0] = rand_dir[0], dir[1] = rand_dir[1], dir[2] = rand_dir[2];

  double pos[3] = {start[0], start[1], start[2]};

  // where am i
  moab::EntityHandle eh = where_am_i(DAG, RayHistory, start, dir);

  // transport loop
  moab::EntityHandle next_surf = 1;
  while (next_surf) {
    // fire ray
    double dist = 0.;
    moab::EntityHandle surf = fire_ray(DAG, RayHistory, eh, start, dir, dist);
    // update particle position
    pos[0] += dir[0] * dist;
    pos[1] += dir[1] * dist;
    pos[2] += dir[2] * dist;
    // update the next surface
    if (!surf)
      return next_surf; // if next surf is 0 we are done
    // used instead of graveyard
    next_surf = surf;
    // update the volume that we are in
    moab::EntityHandle next_vol = next_volume(DAG, surf, eh);

    // reuse the loop variable
    eh = next_vol;
  }
  // done
  return 0;
}

TEST_F(DagmcThreadTest, dagmc_point_in) {
  int nps = 5000;
  double start_pos[3] = {0., 0., 0.};
  omp_set_num_threads(num_threads);
  
  std::map<moab::EntityHandle,int> score;
  // for each surface in the problem initialise the score to 0;
  int num_surf = DTM->get_dagmc_instance(0)->num_entities(2);
  for ( int j = 1 ; j <= num_surf ; j++ ) {
    moab::EntityHandle eh = DTM->get_dagmc_instance(0)->entity_by_index(2, j);
    score[eh] = 0;
  }
  
#pragma omp parallel for shared(score)
  //std::cout << omp_get_num_threads() << std::endl;
  for (int i = 1 ; i <= nps ; i++) {
    //std::cout << i << std::endl;
    moab::EntityHandle eh = transport_cycle(DTM->get_dagmc_instance(omp_get_thread_num()),
					    &(DTM->get_dagmc_raystate(omp_get_thread_num())->history),i, start_pos);
    # pragma omp atomic
    score[eh] +=1;
    // end of transport loop reset state
    DTM->get_dagmc_raystate(omp_get_thread_num())->history.reset();
  }

  int sum = 0;
  for ( int j = 1 ; j <= num_surf ; j++ ) {
    moab::EntityHandle eh = DTM->get_dagmc_instance(0)->entity_by_index(2,j);
    std::cout << eh << " " << score[eh] << std::endl;
    sum += score[eh];
  }
  std::cout << "Total particle leaving " << sum << std::endl;
  std::cout << "Total particles starting " << nps << std::endl;

  delete DTM;
}

TEST_F(DagmcThreadTest, obb_check) {
  // set the number of threads to be the number of procs
  DagThreadManager *dtm = new DagThreadManager(num_threads);
  // get a dagmc instance to continue setup
  moab::DagMC* DAG = dtm->get_dagmc_instance(0);
  moab::ErrorCode rloadval = DAG->load_file(input_file);
  assert(rloadval == moab::MB_SUCCESS);
  // Create the OBB
  rval = DAG->init_OBBTree();
  rval = DAG->moab_instance()->write_mesh("test1.h5m");
  assert(rval == moab::MB_SUCCESS);
  //    DTM->setup_child_threads();
  DTM->initialise_child_threads();
  rval=DAG->moab_instance()->write_mesh("test2.h5m");

  delete dtm;
}
