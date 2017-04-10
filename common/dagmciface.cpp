#include "dagmciface.h"

#include "DagMC.hpp"
using moab::DagMC;


#define DAG DagMC::instance()

void dagmc_version_(double* dagmcVersion) {
  *dagmcVersion = DAG->version();
}

