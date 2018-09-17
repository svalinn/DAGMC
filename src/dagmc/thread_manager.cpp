#include <iostream>
#include "thread_manager.hpp"

// constructor
DagThreadManager::DagThreadManager(int thread_count, moab::Interface* MBI) {
  // if the moab pointer is not null point to it
  if (MBI != NULL) {
    MOAB = MBI;
  } else {
    // make new moab instance to share all DAGMC data
    MOAB = new moab::Core();
  }

  // make a new GTT for Thread manager
  GTT = new moab::GeomTopoTool(MOAB,false);
  
  // number of threads
  set_num_threads(thread_count);

  // set storage for master
  dagmc_instances.push_back(new moab::DagMC(GTT));
  dagmc_rayhistories.push_back(new DagMCRayState());

  // setup threads
  setup_child_threads();
}


DagThreadManager::DagThreadManager(moab::Interface* MBI) {
  // if the moab pointer is not null point to it
  if (MBI != NULL) {
    MOAB = MBI;
  } else {
    // make new moab instance to share all DAGMC data
    MOAB = new moab::Core();
  }

  // make a new GTT
  GTT = new moab::GeomTopoTool(MOAB,false);
  
  // set storage for master
  dagmc_instances.push_back(new moab::DagMC(GTT));
  dagmc_rayhistories.push_back(new DagMCRayState());
}

// destructor
DagThreadManager::~DagThreadManager() {
  for (int i = 0 ; i < num_threads ; i++) {
    // delete each dagmc instance
    delete dagmc_instances[i];
    delete dagmc_rayhistories[i];
  }
}

// set the number of threads
void DagThreadManager::set_num_threads(int thread_count) {
  num_threads = thread_count;
}

void DagThreadManager::setup_child_threads() {
  // loop over the number of threads and make a new DAGMC instance for each
  for (int i = 1 ; i < num_threads ; i++) {
    // push back a vector of DAGMC instances
    dagmc_instances.push_back(new moab::DagMC(GTT));
    // collection of ray histories for each thread
    dagmc_rayhistories.push_back(new DagMCRayState());
  }
}

// initalise threads with the loaded data
void DagThreadManager::initialise_child_threads() {
  for (int i = 1 ; i < num_threads ; i++) {
    moab::ErrorCode rval;
    //std::cout << i << std::endl;
    rval = get_dagmc_instance(i)->load_existing_contents();
    rval = get_dagmc_instance(i)->init_OBBTree();
  }
}
