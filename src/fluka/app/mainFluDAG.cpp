//-------------*-C++, //Fortran-*--------------------//
/*!
 * \file   ~/DAGMC/FluDAG/src/cpp/mainFluDAG.cpp
 * \author Julie Zachman
 * \date   Apr 5 2013
 * \brief  Functions called by fluka
 * \note   Unittested
 */
//----------------------------------------------------//
#include <time.h>  // for timing the routine
#include <fstream>
#include <cstdlib>
#include <cstring>

//---------------------------------------------------------------------------//
#include "fluka_funcs.h"
#include "DagMC.hpp"
#include "dagmcmetadata.hpp"
#include "moab/ProgOptions.hpp"

#define flukam flukam_

moab::DagMC* DAG = new moab::DagMC();  // dagmc instance

#ifdef __cplusplus
extern "C" {
#endif

void flukam(const int& GeoFlag);

#ifdef __cplusplus
}
#endif

// This has been modified so that all runs are now fluka runs.
// Getting the geometry out of the file is done by dagmc_define
int main(int argc, char* argv[]) {
  bool flukarun = false;
  moab::ErrorCode error;
  time_t time_before, time_after;

  // Default h5m filename is for fluka runs
  std::string infile = "dagmc.h5m";
  std::string dagmc_file = "";

  // form the inputs and determine if this is a true calculation or a preprocess run
  ProgOptions po("mainfludag: a DAGMC enabled version of FLUKA-CERN");
  po.addOpt<std::string>("dagmc", "Path to h5m DAGMC file to proccess", &dagmc_file);
  po.addOptionalArgs<std::string>(0,"","");
  po.parseCommandLine(argc, argv);

  // if no string has been provided, dagmc command not applied
  // we assume that its a FLUKA run
  if (dagmc_file.empty()) {
    flukarun = true;
    // fluka creates a run dir one lvl higher than cwd
    infile = "../" + infile;
  } else {
    // dagmc command has been set do the its a process run
    infile = dagmc_file;
  }

  // test to see if the file exists
  std::ifstream h5mfile(infile.c_str());  // filestream for mesh geom
  if (!h5mfile.good()) {
    std::cout << "h5m file does not exist" << std::endl;
    exit(1);
  }

  // get the current time
  time(&time_before);  /* get current time; same as: timer = time(NULL)  */

  // DAG call to load the file
  error = DAG->load_file(infile.c_str()); // load the dag file takeing the faceting from h5m
  if (error != moab::MB_SUCCESS) {
    std::cerr << "DAGMC failed to read input file: " << infile << std::endl;
    exit(EXIT_FAILURE);
  }

  time(&time_after);

  double seconds = difftime(time_after, time_before); //get the time in seconds to load file
  time_before = time_after; // reset time to now for the next call
  std::cout << "Time to load the h5m file = " << seconds << " seconds" << std::endl;

  // DAG call to initialize geometry
  // if more than 1 argument provided
  // this an actual calculation
  if (flukarun) {
    error = DAG->init_OBBTree();
  } else {
    // otherwise this is a preprocess run
    // no need to build the tree - its faster
    error = DAG->setup_impl_compl();
    error = DAG->setup_indices();
  }

  // check
  if (error != moab::MB_SUCCESS) {
    std::cerr << "DAGMC failed to initialize geometry and create OBB tree"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  time(&time_after);
  seconds = difftime(time_after, time_before);
  std::cout << "Time to initialise the geometry" << seconds << std::endl;

  // if fluka preprocess run then create mat file to paste into input deck
  if (!flukarun) {
    std::string lcad = "mat.inp";
    std::cout << "Producing material snippets" << std::endl;
    std::cout << "please paste these into your input deck" << std::endl;
    fludag_write(infile, lcad);

    std::string vol_id = "vol_id_idx.txt";
    std::cout << "Producing volume index & id correspondences" << std::endl;
    fludag_write_ididx(vol_id);
  } else {
    // call fluka run
    
    // check for the input file argument
    // get it from the command line
    if(argc >= 1) {
      // convert to std::string
      std::string chinpf_s(argv[1]);
      char chinpf[256] = "";
      memset(chinpf,' ',256);
      std::copy(chinpf_s.begin(),chinpf_s.end(),chinpf);
      strcpy(chcmpt_.chinpf,chinpf);
    } else {
      // get it from the environment
      std::cout << "from env" << std::endl;
      char* env = std::getenv("INPF");
      std::cout << env << std::endl;
      strncpy(env,chcmpt_.chinpf,sizeof(env));
      if (chcmpt_.chinpf[0] == 0) {
	//	flabrt("FLUKAM","FLUDAG NO INPUT SPECIFIED");
	return 1;
      }
    }
    const int flag = 2;
    std::cout << "running fluka" << std::endl;
    flukam(flag);
  }

  return 0;
}
