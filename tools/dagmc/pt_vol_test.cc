#include "MBInterface.hpp"
#include "MBCore.hpp"
#include "DagMC.hpp"
#include "MBTagConventions.hpp"

#include <vector>
#include <iostream>
#include <math.h>
#include <limits>
#include <stdlib.h>

#define CHKERR if (MB_SUCCESS != rval) return rval



MBErrorCode test_pt_volume(DagMC &dagmc, int volID, double xxx, double yyy, double zzz, int &inside)
{
  MBErrorCode rval;

  MBEntityHandle vol = dagmc.entity_by_id(3,volID);

  double u, v, w;
  rval = dagmc.point_in_volume( vol, xxx, yyy, zzz, inside, u, v, w);
  CHKERR;
  
  return MB_SUCCESS;

}

MBErrorCode test_pt_volume_slow(DagMC &dagmc, int volID, double xxx, double yyy, double zzz, int &inside)
{
  MBErrorCode rval;

  MBEntityHandle vol = dagmc.entity_by_id(3,volID);

  rval = dagmc.point_in_volume_slow( vol, xxx, yyy, zzz, inside);
  CHKERR;
  
  return MB_SUCCESS;

}

int main( int argc, char* argv[] )
{
  MBErrorCode rval;

  if (argc < 6) {
    std::cerr << "Usage: " << argv[0] << " <mesh_filename> "
              << " <vol_id> <xxx> <yyy> <zzz> " << std::endl;
    return 1;
  }
  
  char* filename = argv[1];
  int volID = atoi(argv[2]);
  double xxx = atof(argv[3]);
  double yyy = atof(argv[4]);
  double zzz = atof(argv[5]);
  int inside;
  
  std::cout << "Checking pt_in_volume for:" << std::endl
	    << "(x,y,z) = (" << xxx << "," << yyy << "," << zzz << ")"
	    << std::endl << "in volume " << volID 
	    << " of geometry " << filename << std::endl;

  
  DagMC& dagmc = *DagMC::instance();
  rval = dagmc.load_file_and_init( filename, strlen(filename), 0, 0 );
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to initialize DagMC." << std::endl;
    return 2;
  }
  
  int errors = 0;
  
  rval = test_pt_volume(dagmc,volID,xxx,yyy,zzz,inside);
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to test point in volume [fast]." << std::endl;
    return 3;
  }

  std::cout << "[fast] Point is " << (inside?"inside":"outside") << " volume " 
	    << volID << std::endl;
  
  rval = test_pt_volume_slow(dagmc,volID,xxx,yyy,zzz,inside);
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to test point in volume [slow]." << std::endl;
    return 3;
  }

  std::cout << "[slow] Point is " << (inside?"inside":"outside") << " volume " 
	    << volID << std::endl;
  
  return errors;
}
