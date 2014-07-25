//----------------------------------*-C++, Fortran-*----------------------------------// /*!
/* \file   ~/DAGMC/FluDAG/src/cpp/fluka_funcs.cpp
 * \author Julie Zachman 
 * \date   Mon Mar 22 2013 
 * \brief  Functions called by fluka
 * \note   After mcnp_funcs
 */
//---------------------------------------------------------------------------//
// $Id: 
//---------------------------------------------------------------------------//

#include "fluka_funcs.h"

#include "MBInterface.hpp"
#include "MBCartVect.hpp"

#include "DagMC.hpp"
#include "moab/Types.hpp"

using moab::DagMC;

#include <iomanip>
#include <fstream>     // ofstream
#include <sstream>
#include <cstring>
#include <list>        
#include <algorithm>   // sort
#include <utility>     // makepair
#include <stdlib.h>    // atoi

#ifdef CUBIT_LIBS_PRESENT
#include <fenv.h>
#endif

// globals

#define DAG DagMC::instance()

#define DGFM_SEQ   0
#define DGFM_READ  1
#define DGFM_BCAST 2


#ifdef ENABLE_RAYSTAT_DUMPS


#include <fstream>
#include <numeric>

static std::ostream* raystat_dump = NULL;

#endif 
#define DEBUG 1
/* These 37 strings are predefined FLUKA materials. Any ASSIGNMAt of unique 
 * materials not on this list requires a MATERIAL card. */
std::string flukaMatStrings[] = {"BLCKHOLE", "VACUUM", "HYDROGEN",
"HELIUM", "BERYLLIU", "CARBON", "NITROGEN", "OXYGEN", "MAGNESIU",      
"ALUMINUM", "IRON", "COPPER", "SILVER", "SILICON", "GOLD", "MERCURY",  
"LEAD", "TANTALUM", "SODIUM", "ARGON", "CALCIUM", "TIN", "TUNGSTEN",   
"TITANIUM", "NICKEL", "WATER", "POLYSTYR", "PLASCINT", "PMMA",         
"BONECOMP", "BONECORT", "MUSCLESK", "MUSCLEST", "ADTISSUE", "KAPTON",  
"POLYETHY", "AIR"};

int NUM_FLUKA_MATS = 37;

bool debug = false; 

std::set<int> make_exception_set()
{
    std::set<int> nuc_exceptions;
    
    // Preserve FLUKA Entropy: Ref Fluka Manual pp. 318-323
    // Question about 
    // Xenon (many named isotopes, useful only for detectors?)
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("H-1")));   // HYDROG-1
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("H-2")));   // DEUTERIU
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("H-3")));   // TRITIUM
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("He-3")));  // HELIUM-3
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("He-4")));  // HELIUM-4
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Li-6")));  // LITHIU-6
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Li-7")));  // LITHIU-7
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("B-10")));  // BORON-10
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("B-11")));  // BORON-11
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Sr-90"))); // 90-SR
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("I-129"))); // 129-I
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Cs-135")));// 135-CS
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Cs-137")));// 137-CS
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Th-230")));// 230-TH
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Th-232")));// 232-TH
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("U-233"))); // 233-U
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("U-234"))); // 234-U
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("U-235"))); // 235-U
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("U-238"))); // 238-U
    /** This list consists of elements from flukaMatStrings, but it is
        problematic
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("H")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("He")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Be")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("C")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("N")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("O")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Mg")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Al")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Fe")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Cu")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Ag")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Si")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Au")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Hg")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Pb")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Ta")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Na")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Ar")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Ca")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Sn")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("W")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Ti")));
    nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Ni")));
    */

    // Print out results
    debug = 1;
    if (debug)
    {
       std::cout << "Nucids of FLUKA exceptions" << std::endl;
       for (std::set<int>::iterator ptr = nuc_exceptions.begin(); 
            ptr != nuc_exceptions.end(); ++ptr)
       {
            std::cout << *ptr << ", ";
       }
       std::cout << std::endl;
    }

    return nuc_exceptions;
}

const char *delimiters = ":/";

// an empty synonym map to provide as a default argument to parse_properties()
static const std::map<std::string,std::string> no_synonyms;

/* Create a set out of the hardcoded string array. */
std::set<std::string> FLUKA_mat_set(flukaMatStrings, flukaMatStrings+NUM_FLUKA_MATS); 

/* Maximum character-length of a cubit-named material property */
int MAX_MATERIAL_NAME_SIZE = 32;


/* Static values used by dagmctrack_ */

static DagMC::RayHistory history;

bool on_boundary;
double old_direction[3];
MBEntityHandle next_surf; // the next suface the ray will hit
MBEntityHandle prev_surf; // the last value of next surface
MBEntityHandle PrevRegion; // the integer region that the particle was in previously



/**************************************************************************************************/
/******                                FLUKA stubs                                         ********/
/**************************************************************************************************/
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// jomiwr(..)
//---------------------------------------------------------------------------//
/// Initialization routine, was in WrapInit.c
//  For DAGMC only sets the number of volumes in the problem
void jomiwr(int & nge, const int& lin, const int& lou, int& flukaReg)
{

  if(debug)
    {
      std::cout << "================== JOMIWR =================" << std::endl;
    }

  //Original comment:  returns number of volumes
  unsigned int numVol = DAG->num_entities(3);
  flukaReg = numVol;

  if(debug)
    {
      std::cout << "Number of volumes: " << flukaReg << std::endl;
      std::cout << "================== Out of JOMIWR =================" << std::endl;
    }

  return;
}

//---------------------------------------------------------------------------//
// g_step(..)
//---------------------------------------------------------------------------//
//  returns approved step of particle and all variables 
//
void g_step(double& pSx, 
          double& pSy, 
          double& pSz, 
          double* pV,
          int& oldReg,         // pass through
          const int& oldLttc,  // ignore
          double& propStep,    // .
          int& nascFlag,       // .
          double& retStep,     // reset in this method
          int& newReg,         // return from callee
          double& saf,         // safety 
          int& newLttc,        // .
          int& LttcFlag,       // . 
          double* sLt,         // .
          int* jrLt)           // .
{
  double safety; // safety parameter

  if(debug)
    {
      std::cout<<"============= G_STEP	 =============="<<std::endl;    
      std::cout << "Position " << pSx << " " << pSy << " " << pSz << std::endl;
      std::cout << "Direction vector " << pV[0] << " " << pV[1] << " " << pV[2] << std::endl;
      std::cout << "Oldreg = " << oldReg << std::endl;
      std::cout << "PropStep = " << propStep << std::endl;
    }
  
  double point[3] = {pSx,pSy,pSz};
  double dir[3]   = {pV[0],pV[1],pV[2]};  

  if(debug)
    {
      std::cout << "cel = " << oldReg << " pos = " << point[0] << " " << point[1] << " " << point[2];
      std::cout << " dir = " << dir[0] << " " << dir[1] << " " << dir[2] ;
      std::cout << " prop = " << propStep ;
    }
  g_fire(oldReg, point, dir, propStep, retStep, saf, newReg); // fire a ray 
  old_direction[0]=dir[0],old_direction[1]=dir[1],old_direction[2]=dir[2];
  if(debug)
    {
      std::cout << " ret = " << retStep;
      std::cout << " new cel = " << newReg << std::endl;

      std::cout << "saf = " << saf << std::endl;
      std::cout << std::setw(20) << std::scientific;
      std::cout << "newReg = " << newReg << " retStep = " << retStep << std::endl;
    }

  return;
}

//---------------------------------------------------------------------------//
// void g_fire(int& oldRegion, double point[], double dir[], 
//              double &propStep, double& retStep,  int& newRegion)
//---------------------------------------------------------------------------//
// oldRegion - the region of the particle's current coordinates
// point     - the particle's coordinate location vector
// dir       - the direction vector of the particle's current path (ray)
// propStep  - ??
// retStep   - returned as the distance from the particle's current location, along its ray, to the next boundary
// newRegion - gotten from the value returned by DAG->next_vol
// newRegion is gotten from the volue returned by DAG->next_vol
void g_fire(int &oldRegion, double point[], double dir[], double &propStep, 
            double &retStep, double &safety,  int &newRegion)
{

  MBEntityHandle vol = DAG->entity_by_index(3,oldRegion);
  double next_surf_dist;
  MBEntityHandle newvol = 0;


  /*
  if(!check_vol(point,dir,oldRegion))
    {
      history.reset();
    }
  */
  // direction changed reset history
  if( dir[0] == old_direction[0] && dir[1] == old_direction[1] && dir[2] == old_direction[2] ) 
    {   
    }
  else
    {
      history.reset();
    }



  // 
   
  oldRegion = DAG->index_by_handle(vol); // convert oldRegion int into MBHandle to the volume
  if(on_boundary)
    {
      if(boundary_test(vol,point,dir)==0) // if ray not on leaving vol
	{
	  history.reset(); // reset history
	  on_boundary = false; // reset on boundary
	}
    }
  

  MBErrorCode result = DAG->ray_fire(vol, point, dir, next_surf, next_surf_dist,&history); // fire a ray 
  if ( result != MB_SUCCESS )
    {
      std::cout << "DAG ray fire error" << std::endl;
      exit(0);
    }

  if ( next_surf == 0 ) // if next_surface is 0 then we are lost
    {
      std::cout << "!!! Lost Particle !!! " << std::endl;
      std::cout << "in region, " << oldRegion << " aka " << DAG->entity_by_index(3,oldRegion) << std::endl;  
      std::cout.precision(25);
      std::cout << std::scientific ; 
      std::cout << "position of particle " << point[0] << " " << point[1] << " " << point[2] << std::endl;
      std::cout << " traveling in direction " << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
      std::cout << "!!! Lost Particle !!!" << std::endl;
      newRegion = -3; // return error
      return;
    }

  // set the safety
  
  retStep = next_surf_dist; // the returned step length is the distance to next surf
  if ( propStep >= retStep ) // will cross into next volume next step
    {
      MBErrorCode rval = DAG->next_vol(next_surf,vol,newvol);
      newRegion = DAG->index_by_handle(newvol);
      retStep = retStep; //+1.0e-9 ; // path limited by geometry
      next_surf = next_surf;
      on_boundary=true;
      // history is preserved
    }
  else // step less than the distance to surface
    {
      newRegion = oldRegion; // dont leave the current region
      retStep = propStep; //physics limits step
      next_surf = prev_surf; // still hit the previous surface
      history.reset();
	//_to_last_intersection();
      on_boundary=false;
    }

  PrevRegion = newRegion; // particle will be moving to PrevRegion upon next entry.

  if(debug)
  {
     std::cout << "Region on other side of surface is  = " << newRegion << \
                  ", Distance to next surf is " << retStep << std::endl;
  }

  prev_surf = next_surf; // update the surface

  return;
}
///////			End g_step and g_fire
/////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------//
// normal
//---------------------------------------------------------------------------//
// Local wrapper for fortran-called, f_normal.  This function is supplied for testing
// purposes.  Its signature shows what parameters are being used in our wrapper 
// implementation.  
// Any FluDAG calls to f_normal should use this call instead.
// ASSUMES:  no ray history
// Notes
// - direction is not taken into account 
// - curRegion is not currently used.  
int  normal (double& posx, double& posy, double& posz, double *norml, int& curRegion)
{
   int flagErr; 
   int dummyReg;
   double dummyDirx, dummyDiry, dummyDirz;
   f_normal(posx, posy, posz, dummyDirx, dummyDiry, dummyDirz, norml, curRegion, dummyReg, flagErr);
   return flagErr;
}
//---------------------------------------------------------------------------//
// f_normal(..)
//---------------------------------------------------------------------------//
//  Note:  The normal is calculated at the point on the surface nearest the 
//         given point
// ASSUMES:  Point is on the boundary
// Parameters Set:
//     norml vector
//     flagErr = 0 if ok, !=0 otherwise
// Does NOT set any region, point or direction vector.
// Globals used:
//     next_surf, set by ray_fire 
void f_normal(double& pSx, double& pSy, double& pSz,
            double& pVx, double& pVy, double& pVz,
	    double* norml, const int& oldRegion, 
	    const int& newReg, int& flagErr)
{
  if(debug)
  {
      std::cout << "============ NRMLWR =============" << std::endl;
  }

  MBEntityHandle OldReg = DAG -> entity_by_index(3,oldRegion); // entity handle
  double xyz[3] = {pSx,pSy,pSz}; //position vector
  double uvw[3] = {pVx,pVy,pVz}; //particl directoin
  int result; // particle is entering or leaving

  MBErrorCode ErrorCode = DAG->test_volume_boundary( OldReg, next_surf,xyz,uvw, result, &history);  // see if we are on boundary
  ErrorCode = DAG->get_angle(next_surf,xyz,norml); 
  // result = 1 entering, 0 leaving
  if ( result == 0 ) // vector should point towards OldReg
    {
      norml[0] = norml[0]*-1.0;
      norml[1] = norml[1]*-1.0;
      norml[2] = norml[2]*-1.0;
    }


  if(debug)
  {
      std::cout << "Normal: " << norml[0] << ", " << norml[1] << ", " << norml[2]  << std::endl;
  }
  return;
}
///////			End normal() and f_normal()

/////////////////////////////////////////////////////////////////////
//
// check_vol(..)
//
// Returns either true or false, if the point pos belongs to the region oldRegion
//
//
inline bool check_vol( double pos[3], double dir[3], int oldRegion)
{
  int is_inside; // in volume or not
  // convert region id into entityhandle
  MBEntityHandle volume = DAG->entity_by_index(3, oldRegion); // get the volume by index
  MBErrorCode code = DAG->point_in_volume(volume, pos, is_inside,dir);
  if ( code != MB_SUCCESS)
    {
      std::cout << "Failed in DAG call to get point_in_volume" << std::endl;
    }

  if ( is_inside == 1 ) // we are inside the cell tested
    return true;
  else
    return false;
}


/////////////////////////////////////////////////////////////////////
//---------------------------------------------------------------------------//
// look(..)
//---------------------------------------------------------------------------//
// Testable local wrapper for fortran-called, f_look
// This function signature shows what parameters are being used in our wrapper implementation
// oldRegion is looked at if we are no a boundary, but it is not set.
// ASSUMES:  position is not on a boundary
// RETURNS: nextRegion, the region the given point is in 
int look( double& posx, double& posy, double& posz, double* dir, int& oldRegion)
{
   int flagErr;
   int lattice_dummy;  // not used
   int nextRegion;     
   f_look(posx, posy, posz, dir, oldRegion, lattice_dummy, nextRegion, flagErr, lattice_dummy);
   return nextRegion;
}
//---------------------------------------------------------------------------//
// f_look(..)
//---------------------------------------------------------------------------//
// Wrapper for localisation of starting point of particle.
//
// Question:  Should pV, the direction vector, be used?  
//////////////////////////////////////////////////////////////////
// This function answers the question What volume is the point in?  
// oldReg - not looked at UNLESS the volume is on the boundary, then newReg=oldReg
// nextRegion - set to the volume index the point is in.
// ToDo:  Is there an error condition for the flagErr that is guaranteed not to be equal to the next region?
//        Find a way to make use of the error return from point_in_volume
void f_look(double& pSx, double& pSy, double& pSz,
          double* pV, const int& oldReg, const int& oldLttc,
          int& nextRegion, int& flagErr, int& newLttc)
{
  if(debug)
  {
      std::cout << "======= LKWR =======" << std::endl;
      std::cout << "position is " << pSx << " " << pSy << " " << pSz << std::endl; 
  }
  
  history.reset();

  double xyz[] = {pSx, pSy, pSz};       // location of the particle (xyz)
  const double dir[] = {pV[0],pV[1],pV[2]};
  // Initialize to outside boundary.  This value can be 0 or +/-1 for ouside, inside, or on boundary.
  // ToDo:  Should this be initialized at all?  Or should it be initialized to an invalide number?
  int is_inside = 0;                    
  int num_vols = DAG->num_entities(3);  // number of volumes

  for (int i = 1 ; i <= num_vols ; i++) // loop over all volumes
    {
      MBEntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
      // No ray history 
      MBErrorCode code = DAG->point_in_volume(volume, xyz, is_inside);

      // check for non error
      if(MB_SUCCESS != code) 
	{
	  std::cout << "Error return from point_in_volume!" << std::endl;
	  flagErr = -3;
	  return;
	}

      if ( is_inside == 1 ) // we are inside the cell tested
	{
	  nextRegion = i;
          //BIZARRELY - WHEN WE ARE INSIDE A VOLUME, BOTH, nextRegion has to equal flagErr
	  flagErr = nextRegion;
	  return;	  
	}
      else if ( is_inside == -1 )
	{
	  std::cout << "We cannot be here" << std::endl;
	  exit(0);
	}
    }  // end loop over all volumes

  // if we are here do slow check
  // slow_check(xyz,dir,nextRegion);
  flagErr = nextRegion; // return nextRegion
  return;
}

void f_lostlook(double& pSx, double& pSy, double& pSz,
          double* pV, const int& oldReg, const int& oldLttc,
          int& nextRegion, int& flagErr, int& newLttc)
{
    f_look(pSx,pSy,pSz,pV,oldReg,oldLttc,nextRegion,flagErr,newLttc);
    return;
}

/*
 * entering or leaving, if particle on boundary 
 */
int boundary_test(MBEntityHandle vol, double xyz[3], double uvw[3])
{
  int result;
  MBErrorCode ErrorCode = DAG->test_volume_boundary(vol,next_surf,xyz,uvw, result,&history);  // see if we are on boundary
  return result;
}
//---------------------------------------------------------------------------//
// slow_check(..)
// Not CALLED
//---------------------------------------------------------------------------//
// Helper function
void slow_check(double pos[3], const double dir[3], int &oldReg)
{
  std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
  std::cout << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
  int num_vols = DAG->num_entities(3);  // number of volumes
  int is_inside = 0;
  for (int i = 1 ; i <= num_vols ; i++) // loop over all volumes
    {
      MBEntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
      MBErrorCode code = DAG->point_in_volume(volume, pos, is_inside,dir); 
      if ( code != MB_SUCCESS)
	{
	 std::cout << "Failure from point in volume" << std::endl;
	 exit(0);
	}

      if ( is_inside == 1) // if in volume
	{
	  oldReg = DAG->index_by_handle(volume); //set oldReg
	  std::cout << pos[0] << " " << pos[1] << " " << pos[2] << " " << oldReg << std::endl;
	  return;
	}
    }

  std::cout << "FAILED SLOW CHECK" << std::endl;
  exit(0);
}

/*
 *   Particle localisation when magnetic field tracking is on
 */
void lkmgwr(double& pSx, double& pSy, double& pSz,
            double* pV, const int& oldReg, const int& oldLttc,
	    int& flagErr, int& newReg, int& newLttc)
{
  

    const double xyz[] = {pSx, pSy, pSz}; // location of the particle (xyz)
    int is_inside = 0; // logical inside or outside of volume
    int num_vols = DAG->num_entities(3); // number of volumes

    for (int i = 1 ; i <= num_vols ; i++) // loop over all volumes
      {
	MBEntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
	// No ray history or ray direction.
	MBErrorCode code = DAG->point_in_volume(volume, xyz, is_inside);

	// check for non error
	if(MB_SUCCESS != code) 
	  {
	    std::cout << "Error return from point_in_volume!" << std::endl;
	    flagErr = 1;
	    return;
	  }

	if ( is_inside == 1 ) // we are inside the cell tested
	  {
	    newReg = i;
	    flagErr = i+1;
	    if(debug)
	      {
		std::cout << "point is in region = " << newReg << std::endl;
	      }
	    return;
	  }
      }  // end loop over all volumes

    std::cout << "particle is nowhere!" << std::endl;
    newReg = -100;
    std::cout << "point is not in any volume" << std::endl;
    return;
}

void f_lookdb(double& pSx, double& pSy, double& pSz,
	    double* pV, const int& oldReg, const int& oldLttc,
	    int& newReg, int& flagErr, int& newLttc)
{
  if(debug)
    {
      std::cout<<"============= F_LooKDB =============="<< std::endl;
    }
  //return region number and dummy variables
  newReg=0;   
  newLttc=0;
  flagErr=-1; 

  return;
}


/*
 * f_g1rt
 */
void f_g1rt(void)
{
  if(debug)
    {
      std::cout<<"============ F_G1RT ============="<<std::endl;
    }
    return;
}

// Set DNEAR option if needed
int f_idnr(const int & nreg, const int & mlat) 

{
	

// returns 0 if user doesn't want Fluka to use DNEAR to compute the 
// step (the same effect is obtained with the GLOBAL (WHAT(3)=-1)
// card in fluka input), returns 1 if user wants Fluka always to use DNEAR.

	return 0;
}

///////////////////////////////////////////////////////////////////
// from WrapReg2Name.cc 
//
// Wrapper for getting region name corresponding to given region number
///////////////////////////////////////////////////////////////////
void rg2nwr(const int& mreg, const char* Vname)
{
  std::cerr << "============= RG2NWR ==============" << std::endl;    
  std::cerr << "mreg=" << mreg << std::endl;
  char * vvname;
  region2name(mreg, vvname);
  Vname = vvname;
  std::cerr << "reg2nmwr: Vname " << Vname<< std::endl;  
  return;
}

///////////////////////////////////////////////////////////////////
// from WrapReg.hh 
//
// Wrapper for scoring hits: previous step end-point is taken from 
// history (and compared with fluka region index, flukaReg),
// then the wrapper returns all the information regarding the 
// volume tree, i.e. returns indMother[] array with all the 
// mother volumes index and repMother[] array with all the 
// mother volumes repetition number.   
///////////////////////////////////////////////////////////////////
void rgrpwr(const int& flukaReg, const int& ptrLttc, int& g4Reg,
            int* indMother, int* repMother, int& depthFluka)
{
  std::cerr << "============= RGRPWR ==============" << std::endl;    
  std::cerr << "ptrLttc=" << ptrLttc << std::endl;
  return;
}

///////////////////////////////////////////////////////////////////
// from WrapMag.hh
//
// Wrapper for geometry tracking in magnetic field: returns magnetic 
// field values in a given position.
/////////////////////////////////////////////////////////////////
void fldwr(const double& pX, const double& pY, const double& pZ,
            double& cosBx, double& cosBy, double& cosBz, 
            double& Bmag, int& reg, int& idiscflag)

{
  std::cerr<<"================== MAGFLD ================="<<std::endl;
  return;
}

///////////////////////////////////////////////////////////////////
// from WrapFlgfwr.cc
//
// Wrapper for setting of fluka geometry flag
//////////////////////////////////////////////////////////////////
void flgfwr ( int& flkflg )
{
  std::cerr << "=======FLGFWR =======" << std::endl;
  return;
}

///////////////////////////////////////////////////////////////////
// from WrapLookFX.hh
//
// Wrapper for localisation of particle to fix particular conditions.
// At the moment is the same as WrapLookZ.hh. 
//////////////////////////////////////////////////////////////////
void lkfxwr(double& pSx, double& pSy, double& pSz,
            double* pV, const int& oldReg, const int& oldLttc,
	    int& newReg, int& flagErr, int& newLttc)
{
  std::cerr << "======= LKFXWR =======" << std::endl;

  return;
}

/**************************************************************************************************/
/******                                End of FLUKA stubs                                  ********/
/**************************************************************************************************/

/**
 * Helper function for parsing DagMC properties that are integers.
 * Returns true on success, false if property does not exist on the volume,
 * in which case the result is unmodified.
 * If DagMC throws an error, calls exit().
 */
static bool get_int_prop( MBEntityHandle vol, int cell_id, const std::string& property, int& result ){

  MBErrorCode rval;
  if( DAG->has_prop( vol, property ) ){
    std::string propval;
    rval = DAG->prop_value( vol, property, propval );
    if( MB_SUCCESS != rval ){ 
      std::cerr << "DagMC failed to get expected property " << property << " on cell " << cell_id << std::endl;
      std::cerr << "Error code: " << rval << std::endl;
      exit( EXIT_FAILURE );
    }
    const char* valst = propval.c_str();
    char* valend;
    result = strtol( valst, &valend, 10 );
    if( valend[0] != '\0' ){
      // strtol did not consume whole string
      std::cerr << "DagMC: trouble parsing '" << property <<"' value (" << propval << ") for cell " << cell_id << std::endl;
      std::cerr << "       the parsed value is " << result << ", using that." << std::endl;
    }
    return true;
  }
  else return false;

}

/**
 * Helper function for parsing DagMC properties that are doubles.
 * Returns true on success, false if property does not exist on the volume,
 * in which case the result is unmodified.
 * If DagMC throws an error, calls exit().
 */
static bool get_real_prop( MBEntityHandle vol, int cell_id, const std::string& property, double& result ){

  MBErrorCode rval;
  if( DAG->has_prop( vol, property ) ){
    std::string propval;
    rval = DAG->prop_value( vol, property, propval );
    if( MB_SUCCESS != rval ){ 
      std::cerr << "DagMC failed to get expected property " << property << " on cell " << cell_id << std::endl;
      std::cerr << "Error code: " << rval << std::endl;
      exit( EXIT_FAILURE );
    }
    const char* valst = propval.c_str();
    char* valend;
    result = strtod( valst, &valend );
    if( valend[0] != '\0' ){
      // strtod did not consume whole string
      std::cerr << "DagMC: trouble parsing '" << property <<"' value (" << propval << ") for cell " << cell_id << std::endl;
      std::cerr << "       the parsed value is " << result << ", using that." << std::endl;
    }
    return true;
  }
  else return false;

}


//---------------------------------------------------------------------------//
// FluDAG Material Card  Functions
//---------------------------------------------------------------------------//
// 
//---------------------------------------------------------------------------//
void fludag_write(std::string matfile, std::string lfname)
{
  std::map<int, std::pair< int, pyne::Material> > pyne_map;
  int num_vols = fludag_setup();
  if (num_vols > 0)
  {
     pyne_map = pyne_get_materials(matfile, num_vols);
  }

  // ASSIGNMA records:  returns a list of materials for which COMPOUND cards are needed
  std::ostringstream astr;
  std::list<std::pair<int, pyne::Material> > notBuiltIn;
  notBuiltIn = fludagwrite_assignma(astr, num_vols, pyne_map);
  int num_mats = notBuiltIn.size();
  debug = 1;
  if (debug)
  {
     std::cout << astr.str();
     std::cout << "Number of materials not built-in = " << num_mats << std::endl;
  }

  // Process the matNamesList list so that it is unique
  // matNamesList.sort();
  // matNamesList.unique();

  // Print the final list of unique materials
  /*

  if (debug)
  {
     std::cout << "num_mats  " << num_mats << std::endl;
     std::list<std::string>::iterator it; 
     for (it=matNamesList.begin(); it!=matNamesList.end(); ++it)
     {
        std::cout << *it << std::endl;
     }
  }
  */
  ////////////////////////////////////////////////////////////////////////
  // MATERIAL Cards
  std::cerr << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;

  std::set<int> exception_set = make_exception_set();
  
  std::vector<pyne::Material> compounds;
  std::ostringstream mstr;
  int last_id;
  if (num_mats > 0)
  {
     // Do all the materials at once; mstr will be multiline
     // compounds = fludag_write_material(mstr, last_id, num_mats, matfile);
     compounds = fludag_write_material(mstr, last_id, num_mats, exception_set, notBuiltIn);
  }

  std::cerr << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;

  ////////////////////////////////////////////////////////////////////////
  // COMPOUND Cards
  std::ostringstream cstr;
  std::vector<pyne::Material>::iterator mat_ptr;

   // create a file-reading object to read the fluka names
   std::ifstream fin;
   fin.open("../fluka/doc/el.txt");
   if (!fin.good())
   {
     // exit if file not found
     std::cout << "el.txt should be in ../fluka/doc" << std::endl;
     exit(EXIT_FAILURE);
   }
  // This is needed for the *components* of a COMPOUND material
  // The FLUKa names are not stored
  fluka_name_map = readElements(fin);

  for (mat_ptr = compounds.begin(); mat_ptr != compounds.end(); mat_ptr++)
  {
     fludag_write_compound(cstr, last_id, *mat_ptr);
  }

  ////////////////////////////////////////////////////////////////////////
  // Write all the streams to the input file
  std::string header = "*...+....1....+....2....+....3....+....4....+....5....+....6....+....7...";
  std::ofstream lcadfile (lfname.c_str());
  lcadfile << header << std::endl;
  lcadfile << astr.str();
  lcadfile << header << std::endl;
  lcadfile << mstr.str();
  lcadfile << cstr.str();
  lcadfile.close();
//   std::cout << "Writing lcad file = " << lfname << std::endl; 
// Before opening file for writing, check for an existing file;867

/*
  if( lfname != "lcad" ){
    // Do not overwrite a lcad file if it already exists, except if it has the default name "lcad"
    if( access( lfname.c_str(), R_OK ) == 0 ){
      std::cout << "DagMC: reading from existing lcad file " << lfname << std::endl;
      return; 
    }
  }
*/
}
//---------------------------------------------------------------------------//
// fludag_setup()
//---------------------------------------------------------------------------//
// Get the number of volumes (for looping) and parse the properties in 
// the geometry
int fludag_setup()
{
  int num_vols = DAG->num_entities(3);

  MBErrorCode ret;

  std::vector< std::string > keywords;
  ret = DAG->detect_available_props( keywords, delimiters );
  // parse data from geometry so that property can be found
  ret = DAG->parse_properties( keywords, no_synonyms, delimiters );
  if (MB_SUCCESS != ret) 
  {
    std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }
  // Preprocessing loop:  make a string, "props",  for each vol entity
  // This loop could be removed if the user doesn't care to see terminal output
  std::cout << "Property list: " << std::endl;
  for (unsigned i=1; i<=num_vols; i++)
  {
      std::string props = mat_property_string(i, keywords);
      int id = DAG->id_by_index(3, i);
      if (props.size()) 
      {
         std::cout << "Vol " << i << ", id=" << id << ": parsed props: " << props << std::endl; 
      }
      else
      {
         std::cout << "Vol " << i << ", id=" << id << " has no props: " <<  std::endl; 
      }
  }
  return num_vols;
}
//---------------------------------------------------------------------------//
// fludagwrite_assignma
//---------------------------------------------------------------------------//
// Put the ASSIGNMAt statements in the output ostringstream
// std::list<std::string> fludagwrite_assignma(std::ostringstream& ostr, int num_vols, 
// ToDo: Verify this should be a list and not a vector

std::list<std::pair<int, pyne::Material> >fludagwrite_assignma(std::ostringstream& ostr, int num_vols, 
                                          std::map<int, std::pair< int, pyne::Material> > pyne_map)  
{
  std::list<std::pair<int, pyne::Material> > notBuiltIn;

  std::map<int, std::pair< int, pyne::Material> >::iterator mptr; 
  for (mptr=pyne_map.begin(); mptr != pyne_map.end(); ++mptr)
  {
      int id = mptr->first;  
      std::pair<int, pyne::Material> fluka_mat_id_pair = mptr->second;
      pyne::Material mat = fluka_mat_id_pair.second;
     
      // pyne::Material mat = mptr->second.second;
      std::string fluka_name = mat.metadata["fluka_name"].asString();
      if (FLUKA_mat_set.find(fluka_name) == FLUKA_mat_set.end())
      { 
          // current material is not in the pre-existing FLUKA material list
          notBuiltIn.push_back(fluka_mat_id_pair); 
      }
      ostr << std::setw(10) << std::left  << "ASSIGNMAt";
      ostr << std::setw(10) << std::right << fluka_name;
      ostr << std::setprecision(0) << std::fixed << std::showpoint 
         << std::setw(10) << std::right << (float)id << std::endl;
      
  }
//   std::map<int, int> id_index_map = write_id_index(num_vols);

  // Finish the ostr with the implicit complement card
  std::string implicit_comp_comment = "* The next volume is the implicit complement";
  ostr << implicit_comp_comment << std::endl;
  ostr << std::setw(10) << std::left  << "ASSIGNMAt";
  ostr << std::setw(10) << std::right << "VACUUM";
  ostr << std::setprecision(0) << std::fixed << std::showpoint 
         << std::setw(10) << std::right << (float)num_vols << std::endl;
  debug = 1;
  if (debug)
  {
     // Show the output string just created
     std::cout << ostr.str();
  }
  ////////////////////////////////////////////////////////////////////////////////
  // 
  // At this point we have 
  //  o all the ASSIGNMa cards collected in ostr
  //  o a list of material names (that can be made unique) 
  //    to be used for the material cards
  //  o written a helpuful file out called "index_id.txt" 
  ////////////////////////////////////////////////////////////////////////////////

  return notBuiltIn;
}

/*
std::map<int, int> write_id_index(int num_vols)
{
  std::map<int, int> id_index_map;
  // Open an outputstring for index-id table and put a header in it
  std::ostringstream idstr;
  idstr << std::setw(5) <<  "Index" ;
  idstr << std::setw(5) <<  "   Id" << std::endl;
  MBEntityHandle entity = 0;
  // Create the index_id map
  for (unsigned vol_i = 0 ; vol_i <= num_vols ; vol_i++)
  {
      entity = DAG->entity_by_index(3, vol_i);
      int id = DAG->id_by_index(3, vol_i);

      idstr << std::setw(5) << std::right << vol_i;
      // Mark the vol index negative if it's the graveyard
      // ToDo:  This loop is never hit.  There is a material named "Graveyar"
      if (DAG->has_prop(entity, "graveyard"))
      {
         std::cerr << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
         id = -1;
      }
      idstr << std::setw(5) << std::right << id << std::endl;
      id_index_map.insert(std::make_pair(vol_i,id));
  }
  // Prepare an output file for idstr
  std::string index_id_filename = "index_id.txt";
  std::ofstream index_id(index_id_filename.c_str());
  index_id << idstr.str();
  index_id.close(); 

  return id_index_map;
}
*/
//---------------------------------------------------------------------------//
// pyne_get_materials
//---------------------------------------------------------------------------//
// Return the set of all PyNE material objects in the current geometry
// This function mimics what MaterialLibary.from_hdf5(..) might return
std::map<int, std::pair<int, pyne::Material> > pyne_get_materials(std::string mat_file, int num_mats)
{
  pyne::Material mat;
  std::map<int, std::pair< int, pyne::Material> > id_mat_map;

  // Skip the implicit complement
  // Don't look at the graveyard or the material with no properties
  // ToDo:  may want to get a better looping number a priori
  std::cerr << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  for (int i=0; i<num_mats-1; i++)
  {
      mat = pyne::Material();
      mat.from_hdf5(mat_file, "/materials", i, 1);

      if (mat.metadata.isMember("fluka_name"))
      {
         std::string fluka_id_str = mat.metadata["fluka_material_index"].asString();
         std::string mat_id_str = mat.metadata["mat_number"].asString();

         int id = atoi(mat_id_str.c_str());
         int fluka_mat_id = atoi(fluka_id_str.c_str());
         std::cout << id << ". " << mat.metadata["fluka_name"].asString() 
                << " --> " << fluka_mat_id << std::endl;

         id_mat_map[id] = std::make_pair(fluka_mat_id,mat);
      }
      else
      {
         std::cout << "This material does not have fluka_name as a member" << std::endl;
	 print_material(mat);
      }
      /*
      if (mat.metadata.isMember("fluka_name"))
      {
         std::cout << mat.metadata["fluka_name"].asString();
      }
      if (mat.metadata.isMember("comment"))
      {
         std::cout << mat.metadata["comment"].asString();
         std::cout << mat.metadata["comment"];
      }
      
      else
      {
         std::cout << "\t----- no comments ----- ";
      }
     
      std::cout << std::endl; 
      */
  }
  std::cerr << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  return id_mat_map;
}

//---------------------------------------------------------------------------//
// fludag_write_material
//---------------------------------------------------------------------------//
// Return a FLUKA MATERIAL card based on the info that has previously been
// captured from the material file.
// 
std::vector<pyne::Material> fludag_write_material(std::ostringstream& ostr, 
                                                  int& last_id,
                                                  int num_mats, 
                                                  std::set<int> exception_set,
                                                  std::list<std::pair<int, pyne::Material> > notBuiltIn)
{
  pyne::Material mat;
  pyne::Material collmat;
  // For sorting the MATERIAL cards in order of fluka_mat_id
  std::list<std::pair<int, std::string> > id_line_list;
  std::vector<pyne::Material> materials;

  std::list<std::pair<int, pyne::Material> >::iterator ptr;
  for (ptr=notBuiltIn.begin(); ptr != notBuiltIn.end(); ++ptr)
  {
      pyne::Material mat = ptr->second;
      collmat = mat.collapse_elements(exception_set);

      collmat.density = mat.density;
      collmat.metadata = mat.metadata;
      if (debug)
      {
         print_material(collmat);
      }
      collmat = mat.collapse_elements(exception_set);
      collmat.density = mat.density;
      collmat.metadata = mat.metadata;
      if (debug)
      {
         print_material(collmat);
      }
      // Construct the material card.
      std::string line = collmat.fluka();

      // Extract the material id so the lines can be sorted on it
      std::string id_str = mat.metadata["fluka_material_index"].asString();
      int id = atoi(id_str.c_str());
      id_line_list.push_back(std::make_pair(id, line));

      // To prepare for for a potential compound card, add this
      // material to a container if it's not a FLUKA material
      std::string mat_fluka_name = collmat.metadata["fluka_name"].asString();
      if (FLUKA_mat_set.find(mat_fluka_name) == FLUKA_mat_set.end())
      {
         // current material is not in the pre-existing FLUKA material list
         materials.push_back(collmat);
      }
  }

  // Lists magically know how to sort <int, string> pairs
  id_line_list.sort();
  last_id = id_line_list.back().first;

  // Now that its sorted grab the lines sequentially
  while ( 0 < id_line_list.size() )
  {
      ostr << id_line_list.front().second;
      id_line_list.pop_front();
  }
  return materials;
}
//---------------------------------------------------------------------------//
// fludag_write_compound
//---------------------------------------------------------------------------//
// Only material objects that have non-fluka names are processed
struct CompoundElement 
{
   int za_id;              // znum+anum 
   std::string fluka_name; // may or may not be predefined
   bool predefined;        // True if fluka_name is predefined
   std::string name;
   double frac;
};

void fludag_write_compound(std::ostringstream& cstr, int& last_id, 
                           pyne::Material& material)
{
   // create a file-reading object to read the fluka names
   std::ifstream fin;
   fin.open("../fluka/src/el.txt");
   if (!fin.good())
   {
     // exit if file not found
     std::cout << "el.txt should be in ../fluka/src" << std::endl;
     exit(EXIT_FAILURE);
   }

   std::cerr << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
   std::list< CompoundElement > list_ce;

   std::string material_name = material.metadata["fluka_name"].asString();
   std::cout << "Writing compound card for material " << 
                 material.metadata["fluka_name"].asString() << std::endl;

   std::cout << "  Nucid,   Frac,  PyNE name " << std::endl;

   std::map<int, double>::iterator comp_iter = material.comp.begin(); 
   // Go through the composition of this material
   for (; comp_iter != material.comp.end(); comp_iter++)
   {
       CompoundElement ce;
       std::cout << comp_iter->first << ", " << comp_iter->second << ", ";
       std::cout << pyne::nucname::name(comp_iter->first) << std::endl;
       
       int znum = pyne::nucname::znum(comp_iter->first);
       int anum = pyne::nucname::anum(comp_iter->first);
       // Make the za_id
       ce.za_id = znum*10000000 + anum*10000;

       // Use za_id to get the fluka_name from special database
       std::string cur_name = fluka_name_map[ce.za_id];

       ce.fluka_name = cur_name; 
       ce.frac = comp_iter->second;
       if (FLUKA_mat_set.find(cur_name) == FLUKA_mat_set.end())
       {
          // current name is not in the pre-existing FLUKA material list
	  // modify it for the compound card
          // ce.name = modify_fluka_name(cur_name);
	  // jcz test for segmentation
          ce.name = cur_name;
          ce.predefined = false;
       }
       else
       {
          // It is predefined, copy it into the right spot
          ce.name = cur_name;
          ce.predefined = true;
       }
       std::cout << std::endl;

       list_ce.push_back(ce);

   }  // end compounds

   std::string header = "*...+....1....+....2....+....3....+....4....+....5....+....6....+....7...";
   std::cout << header << std::endl;
   // Note: List is popped 3 times each time in loop
   while (list_ce.size() >= 3)
   {
      CompoundElement ce1 = list_ce.front();
      if (!ce1.predefined)
      {
         define_material(cstr, ++last_id, ce1.name);
      }
      list_ce.pop_front();

      CompoundElement ce2 = list_ce.front();
      if (!ce2.predefined)
      {
         define_material(cstr, ++last_id, ce2.name);
      }
      list_ce.pop_front();

      CompoundElement ce3 = list_ce.front();
      if (!ce3.predefined)
      {
         define_material(cstr, ++last_id, ce3.name);
      }
      list_ce.pop_front();

      cstr << std::setw(10) << std::left  << "COMPOUND";
      cstr << std::setw(10) << std::right << ce1.frac;
      cstr << std::setw(10) << std::right << ce1.name;
      cstr << std::setw(10) << std::right << ce2.frac;
      cstr << std::setw(10) << std::right << ce2.name;
      cstr << std::setw(10) << std::right << ce3.frac;
      cstr << std::setw(10) << std::right << ce3.name;
      cstr << std::setw(10) << std::left << material_name;
      cstr << std::endl;
   }
   // Get the last one or two fractions
   if (list_ce.size() > 0)
   {
      CompoundElement ce = list_ce.front();
      if (!ce.predefined)
      {
         define_material(cstr, ++last_id, ce.name);
      }
      list_ce.pop_front();

      cstr << std::setw(10) << std::left  << "COMPOUND";
      cstr << std::setw(10) << std::right << ce.frac;
      cstr << std::setw(10) << std::right << ce.name;
     
      if (list_ce.size() > 0)
      {
         CompoundElement ce = list_ce.front();
         if (!ce.predefined)
         {
            define_material(cstr, ++last_id, ce.name);
         }
         cstr << std::setw(10) << std::right << ce.frac;
         cstr << std::setw(10) << std::right << ce.name;
      }
      else
      {
	 cstr << std::setw(10) << std::right << ""; 
	 cstr << std::setw(10) << std::right << ""; 
      }
      cstr << std::setw(10) << std::right << ""; 
      cstr << std::setw(10) << std::right << ""; 

      cstr << std::setw(10) << std::left << material_name;
   }
}  // end fludag_write_compound

//---------------------------------------------------------------------------//
// define_material
//---------------------------------------------------------------------------//
// Convenience function to create a material with an irrelevant density
// The material is to be part of a compound
void define_material(std::ostringstream &cstr, int &last_id, std::string dname)
{
    std::cout << "; in define_material, last_id = " << last_id << std::endl;
    std::stringstream id_stream;
    id_stream << last_id << ".";
    // Use the defined name in a MATERIAL card
    cstr << std::setw(10) << std::left << "MATERIAL";
    cstr << std::setw(10) << std::right << "";
    cstr << std::setw(10) << std::right << "";
    cstr << std::setprecision(0) << std::fixed << std::showpoint 
          << std::setw(10) << std::right << 999.;
    cstr << std::setw(10) << std::right << id_stream.str();
    cstr << std::setw(10) << std::right << "";
    cstr << std::setw(10) << std::right << "";
    cstr << std::setw(10) << std::left << dname << std::endl;
}

//---------------------------------------------------------------------------//
// modify_fluka_name
//---------------------------------------------------------------------------//
// jcz look at for segmentation
std::string modify_fluka_name(std::string& origName)
{
   if (6 >= origName.length())
   {
      return origName.append("_1");
   }
   if (7 <= origName.length() && origName.length() <= 8)
   {
      return origName.substr(0,6).insert(6,"_1");
   }
}
//---------------------------------------------------------------------------//
// readElements
//---------------------------------------------------------------------------//
// 

const int MAX_CHARS_PER_LINE = 512;
const int MAX_TOKENS_PER_LINE = 20;
const char* const DELIMITER = " ";

// Hardcoded to read a file of fluka-named elements, atomic #, 
// each with 12 or more id's
std::map<int, std::string> readElements(std::ifstream& fin)
{
   int lastNum = 0;
   std::map<int, std::string> map_fluka_name;
   int i = 0;

   // read each line of the file
   while (!fin.eof())
   {
     // read an entire line into memory
     char buf[MAX_CHARS_PER_LINE];
     fin.getline(buf, MAX_CHARS_PER_LINE);

     // parse the line into blank-delimited tokens
     const char* token;

     // parse the line for the first word
     token = strtok(buf, DELIMITER);
     
     if (token) // line not blank
     {
       std::string name = std::string(token);

       // Get the second token here because we'll need the atomic 
       // # either for assignment or for checking
       unsigned int z_num = atoi(strtok(0,DELIMITER));

       // And get the third value, the id
       int id = atoi(strtok(0,DELIMITER));

       int za_id;
       if (z_num == lastNum ) 
       {
          // We have a named isotope, the id is  important
	  za_id = z_num*10000000 + id*10000;
       }
       else
       {
          // Elemental form
	  za_id = z_num*10000000;
	  std::cout << "I'm here" << std::endl;
       }  // end check for duplicate ids
       std::pair< int, std::string>  couplet = std::make_pair(za_id, name);
       std::cout << i++ << ".  " << couplet.first << ", " << couplet.second << std::endl;
       // map_fluka_name.insert( std::make_pair(za_id, name) );
       map_fluka_name.insert(couplet);
       // jcz this is for catching isotopes and ensuring the anum becomes part of 
       // the nucid; NOTE:  this FAILS for uranium and thorium!! and check others
       // also plutonium, and americium
       lastNum = z_num;

     } // end if line not blank
   }   // end while over lines in file
   
   return map_fluka_name;
}

//---------------------------------------------------------------------------//
// print_material
//---------------------------------------------------------------------------//
// Convenience function
void print_material( pyne::Material material) 
{
  pyne::comp_iter it;

  std::cout << "density = " << material.density << std::endl;
  std::cout << "mass = " << material.mass << std::endl;
  std::cout << "atoms_per_mol = " << material.atoms_per_molecule << std::endl;

  for ( it = material.comp.begin() ; it != material.comp.end() ; ++it ) {
    if(it->second <= 0.0)
      continue;
    else
      std::cout << it->first << " " << it->second << std::endl;
  }
  // std::cout << "Metadata:" << std::endl;
  // std::cout << material.metadata << std::endl;
}

//---------------------------------------------------------------------------//
// mat_property_string
//---------------------------------------------------------------------------//
// For a given volume, find the values of all properties named "MAT".   
// Create a string with these properites
// Modified from make_property_string
// This function helps with debugging, but is not germane to writing cards.
std::string mat_property_string (int index, std::vector<std::string> &properties)
{
  MBErrorCode ret;
  std::string propstring;
  MBEntityHandle entity = DAG->entity_by_index(3,index);
  int id = DAG->id_by_index(3, index);
  for (std::vector<std::string>::iterator p = properties.begin(); p != properties.end(); ++p)
  {
     if ( DAG->has_prop(entity, *p) )
     {
        std::vector<std::string> vals;
        ret = DAG->prop_values(entity, *p, vals);
        if( ret != MB_SUCCESS )
        {
            std::cerr << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
            std::exit( EXIT_FAILURE );
        }
        propstring += *p;
        if (vals.size() == 1)
        {
 	   propstring += "=";
           propstring += vals[0];
        }
        else if (vals.size() > 1)
        {
 	   // this property has multiple values; list within brackets
           propstring += "=[";
	   for (std::vector<std::string>::iterator i = vals.begin(); i != vals.end(); ++i)
           {
	       propstring += *i;
               propstring += ",";
           }
           // replace the last trailing comma with a close bracket
           propstring[ propstring.length() -1 ] = ']';
        }
        propstring += ", ";
     }
  }
  if (propstring.length())
  {
     propstring.resize( propstring.length() - 2); // drop trailing comma
  }
  return propstring;
}

//---------------------------------------------------------------------------//
// make_property_string
//---------------------------------------------------------------------------//
// For a given volume, find all properties associated with it, and any and all 
//     values associated with each property
// Copied and modified from obb_analysis.cpp
static std::string make_property_string (MBEntityHandle eh, std::vector<std::string> &properties)
{
  MBErrorCode ret;
  std::string propstring;
  for (std::vector<std::string>::iterator p = properties.begin(); p != properties.end(); ++p)
  {
     if ( DAG->has_prop(eh, *p) )
     {
        std::vector<std::string> vals;
        ret = DAG->prop_values(eh, *p, vals);
        if( ret != MB_SUCCESS )
        {
            std::cerr << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
            std::exit( EXIT_FAILURE );
        }
        propstring += *p;
        if (vals.size() == 1)
        {
 	   propstring += "=";
           propstring += vals[0];
        }
        else if (vals.size() > 1)
        {
 	   // this property has multiple values; list within brackets
           propstring += "=[";
	   for (std::vector<std::string>::iterator i = vals.begin(); i != vals.end(); ++i)
           {
	       propstring += *i;
               propstring += ",";
           }
           // replace the last trailing comma with a close bracket
           propstring[ propstring.length() -1 ] = ']';
        }
        propstring += ", ";
     }
  }
  if (propstring.length())
  {
     propstring.resize( propstring.length() - 2); // drop trailing comma
  }
  return propstring;
}

//////////////////////////////////////////////////////////////////////////
/////////////
/////////////		region2name - modified from dagmcwrite
/////////////
//                      Called in WrapReg2Name
//////////////////////////////////////////////////////////////////////////
void region2name(int volindex, char *vname )  // file with cell/surface cards
{
  MBErrorCode rval;

  std::vector< std::string > fluka_keywords;
  fluka_keywords.push_back( "mat" );
  fluka_keywords.push_back( "rho" );
  fluka_keywords.push_back( "comp" );
  fluka_keywords.push_back( "graveyard" );

  // parse data from geometry
  rval = DAG->parse_properties (fluka_keywords, no_synonyms, delimiters);
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  int cmat = 0;
  double crho;

  MBEntityHandle vol = DAG->entity_by_index( 3, volindex );
  int cellid = DAG->id_by_index( 3, volindex);

  bool graveyard = DAG->has_prop( vol, "graveyard" );

  std::ostringstream istr;
  if( graveyard )
  {
     istr << "BLCKHOLE";
     if( DAG->has_prop(vol, "comp") )
     {
       // material for the implicit complement has been specified.
       get_int_prop( vol, cellid, "mat", cmat );
       get_real_prop( vol, cellid, "rho", crho );
       std::cout 
            << "Detected material and density specified for implicit complement: " 
	    << cmat <<", " << crho << std::endl;
     }
   }
   else if( DAG->is_implicit_complement(vol) )
   {
      istr << "mat_" << cmat;
      if( cmat != 0 ) istr << "_rho_" << crho;
   }
   else
   {
      int mat = 0;
      get_int_prop( vol, cellid, "mat", mat );

      if( mat == 0 )
      {
        istr << "0";
      }
      else
      {
        double rho = 1.0;
        get_real_prop( vol, cellid, "rho", rho );
        istr << "mat_" << mat << "_rho_" << rho;
      }
   }
   char *cstr = new char[istr.str().length()+1];
   std:strcpy(cstr,istr.str().c_str());
   vname = cstr;
}
// Not Called
void dagmc_version_(double* dagmcVersion)
{
  *dagmcVersion = DAG->version();
}
