#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "DagMC.hpp"
#include "MBTagConventions.hpp"

#include <vector>
#include <iostream>
#include <math.h>
#include <limits>
#include <stdlib.h>

#define CHKERR if (MB_SUCCESS != rval) return rval

using namespace moab;

ErrorCode test_pt_volume(DagMC &dagmc, int volID, double xxx, double yyy, double zzz, int &inside,
			 double uuu, double vvv, double www)
{
  ErrorCode rval;

  EntityHandle vol = dagmc.entity_by_id(3,volID);
  double xyz[3] = {xxx,yyy,zzz};
  double uvw[3] = {uuu,vvv,www};

  rval = dagmc.point_in_volume( vol, xyz, inside, uvw );
  CHKERR;
  
  return MB_SUCCESS;

}

ErrorCode test_pt_volume_slow(DagMC &dagmc, int volID, double xxx, double yyy, double zzz, int &inside)
{
  ErrorCode rval;

  EntityHandle vol = dagmc.entity_by_id(3,volID);
  double xyz[3] = {xxx,yyy,zzz};

  rval = dagmc.point_in_volume_slow( vol, xyz, inside);
  CHKERR;
  
  return MB_SUCCESS;

}

int main( int argc, char* argv[] )
{
  ErrorCode rval;

  if (argc != 6 && argc != 9) {
    std::cerr << "Usage: " << argv[0] << " <mesh_filename> "
              << " <vol_id> <xxx> <yyy> <zzz> [<uuu> <vvv> <www>]" << std::endl;
    return 1;
  }
  
  char* filename = argv[1];
  int volID = atoi(argv[2]);
  double xxx = atof(argv[3]);
  double yyy = atof(argv[4]);
  double zzz = atof(argv[5]);
  double uuu = 0;
  double vvv = 0;
  double www = 0;

  if (argc > 6) {
    uuu = atof(argv[6]);
    vvv = atof(argv[7]);
    www = atof(argv[8]);
  }

  int inside;
  
  std::cout << "Checking pt_in_volume for:" << std::endl
	    << "(x,y,z) = (" << xxx << "," << yyy << "," << zzz << ")"
	    << std::endl << "in volume " << volID 
	    << " of geometry " << filename << std::endl;

  
  DagMC dagmc = DagMC();
  rval = dagmc.load_file( filename, 0 );
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to load file." << std::endl;
    return 2;
  }
  rval = dagmc.init_OBBTree( );
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to initialize DagMC." << std::endl;
    return 2;
  }
  
  int errors = 0;
  
  rval = test_pt_volume(dagmc,volID,xxx,yyy,zzz,inside,uuu,vvv,www);
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
