#include "DagMC.hpp"
#include "moab/Core.hpp"

#ifndef THREAD_MANAGER_HPP
#define THREAD_MANAGER_HPP 1

// Struct for storage of ray state
class DagMCRayState {
 public:
  // set some default values via constructor
  DagMCRayState() {
    last_uvw[0] = 0.;
    last_uvw[1] = 0.;
    last_uvw[2] = 0.;
    visited_surface = false;
    use_dist_limit = false;
    double dist_limit = 0.0;
    last_nps = 0;
  }
  ~DagMCRayState() {};

 public:
  double last_uvw[3];
  moab::DagMC::RayHistory history;
  std::vector<moab::DagMC::RayHistory> history_bank;
  std::vector<moab::DagMC::RayHistory> pblcm_history_stack;
  bool visited_surface;
  bool use_dist_limit;
  double dist_limit;
  int last_nps;

};

class DagThreadManager {
 public:
  DagThreadManager(int num_threads, moab::Interface* MBI = NULL);
  DagThreadManager(moab::Interface* MBI = NULL);
  ~DagThreadManager();

  // setup the dagmc state for the threads
  void setup_child_threads();

  void initialise_child_threads();

  void set_num_threads(int thread_count);

  // get the dagmc instance for a given thread
  inline moab::DagMC* get_dagmc_instance(int thread_id) {
    return dagmc_instances[thread_id];
  }

  // get the dagmc_raystate for a given thread
  inline DagMCRayState* get_dagmc_raystate(int thread_id) {
    return dagmc_rayhistories[thread_id];
  }

 private:
  int num_threads; ///< Number of threads to be held by the manager
  std::vector<moab::DagMC*> dagmc_instances; ///< vector of dagmc instances
  std::vector<DagMCRayState*> dagmc_rayhistories; ///< vector to the associated ray history
  moab::Interface* MOAB; ///< moab pointer
  moab::GeomTopoTool *GTT; ///< GTT Pointer
};

#endif
