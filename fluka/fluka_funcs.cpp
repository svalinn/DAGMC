#include "fluka_funcs.h"

#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"
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

#include <fstream>
#include <numeric>

static std::ostream* raystat_dump = NULL;

#define ID_START 26

bool debug = false;

// delimiters used by uwuw
const char *delimiters = ":/";
// an empty synonym map to provide as a default argument to parse_properties()
static const std::map<std::string,std::string> no_synonyms;

/* Maximum character-length of a cubit-named material property */
int MAX_MATERIAL_NAME_SIZE = 32;

// current state of the particle
static particle_state state;

/* For DAGMC only sets the number of volumes in the problem */
void jomiwr(int & nge, const int& lin, const int& lou, int& flukaReg)
{

  if(debug) {
    std::cout << "================== JOMIWR =================" << std::endl;
  }

  //Original comment:  returns number of volumes
  unsigned int numVol = DAG->num_entities(3);
  flukaReg = numVol;

  if(debug) {
    std::cout << "Number of volumes: " << flukaReg << std::endl;
    std::cout << "================== Out of JOMIWR =================" << std::endl;
  }

  return;
}

/* returns the approved step from pS along pV */
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

  if(debug) {
    std::cout<<"============= G_STEP	 =============="<<std::endl;
    std::cout << "Position " << pSx << " " << pSy << " " << pSz << std::endl;
    std::cout << "Direction vector " << pV[0] << " " << pV[1] << " " << pV[2] << std::endl;
    std::cout << "Oldreg = " << oldReg << std::endl;
    std::cout << "PropStep = " << propStep << std::endl;
  }

  double point[3] = {pSx,pSy,pSz};
  double dir[3]   = {pV[0],pV[1],pV[2]};

  g_fire(oldReg, point, dir, propStep, retStep, saf, newReg); // fire a ray

  if(debug) {
    std::cout << " ret = " << retStep;
    std::cout << " new cel = " << newReg << std::endl;
    std::cout << "saf = " << saf << std::endl;
    std::cout << std::setw(20) << std::scientific;
    std::cout << "newReg = " << newReg << " retStep = " << retStep << std::endl;
  }

  return;
}

/* function to determine the particles next volume & step length */
void g_fire(int &oldRegion, double point[], double dir[], double &propStep,
            double &retStep, double &safety,  int &newRegion)
{
  moab::EntityHandle vol = DAG->entity_by_index(3,oldRegion); // get eh of current region
  moab::EntityHandle next_surf; // next surf we hit
  double next_surf_dist;
  moab::EntityHandle newvol = 0;

  if(debug) print_state(state);

  // direction changed reset history, may not be robust
  if( dir[0] != state.old_direction[0] || dir[1] != state.old_direction[1] || dir[2] != state.old_direction[2] )
    if( point[0] != state.old_position[0] || point[1] != state.old_position[1] || point[2] != state.old_position[2] ) {
      // if point and direction has changed, we should reset everything
      reset_state(state);
    } else {
      // direction has changed reset ray history
      state.history.reset();
    }

  /*
  // direction changed reset history, may not be robust
  if( dir[0] == state.old_direction[0] && dir[1] == state.old_direction[1] && dir[2] == state.old_direction[2] )  {
    // position may have changed since last time, may not be robust
    if( point[0] == state.old_position[0] && point[1] == state.old_position[1] && point[2] == state.old_position[2] ) {
    } else {
      // if point and direction has changed, we should reset everything
      reset_state(state);
    }
  } else {
    // direction has changed reset ray history
    state.history.reset();
  }
  */

  if(debug) print_state(state);

  if(state.on_boundary) {
    // electrons do funny things on the boundary, since electrons can interact on the boundary, it may be the case
    // that the when on the boundary the electron changes direction, not end up in the new region as originally thought
    // when this happens
    // Alfredo says that we only need to check the dot product of the normal of the last surface crossed &
    // the particle direction if less than 0.0 then we return a step of 0.0 and set newRegion to -3

    // since the fluka normal points the other way
    if(state.prev_surface == 0 && state.next_surface != 0) {
      if(debug) std::cout << "Checking to see if we have changed direction using next surface" << std::endl;
      if(dot_product(state.next_surface,point,dir) < 0.0 ) {
        retStep = 0.0;
        newRegion = -3;
        return;
      }
    } else if(state.prev_surface != 0 ) {
      if(debug) std::cout << "Checking to see if we have changed direction using last surface" << std::endl;
      if(dot_product(state.prev_surface,point,dir) < 0.0 ) {
        retStep = 0.0;
        newRegion = -3;
        return;
      }
    } else {
      fludag_abort("g_fire","Failure from within FluDAG, next_surf & prev_surf are both 0",1);
    }

    // if we are not on the boundary but think we should be
    if(boundary_test(vol,point,dir) == 0) { // if ray not on boundary of leaving vol
      if(debug) std::cout << " ray is on boundary" << std::endl;
      state.history.reset(); // reset history
      state.on_boundary = false; // reset on boundary
    }
  }

  // perform the actual ray fire
  moab::ErrorCode rval = DAG->ray_fire(vol, point, dir, next_surf, next_surf_dist, &state.history); // fire a ray
  if ( rval != moab::MB_SUCCESS ) fludag_abort("g_fire","Failure from within DAGMC, ray_fire returned errorcode",rval);

  if ( next_surf == 0 ) { // if next_surface is 0 then we are lost
    if(debug) {
      std::cout << "!!! Lost Particle !!! " << std::endl;
      std::cout << "in region, " << oldRegion << " aka " << DAG->entity_by_index(3,oldRegion) << std::endl;
      std::cout.precision(25);
      std::cout << std::scientific ;
      std::cout << "position of particle " << point[0] << " " << point[1] << " " << point[2] << std::endl;
      std::cout << " traveling in direction " << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
      std::cout << "!!! Lost Particle !!!" << std::endl;
    }
    newRegion = -33; // return error
    retStep = 0.0;
    return;
  }

  // set the safety
  retStep = next_surf_dist; // the returned step length is the distance to next surf

  // proposed step
  double proposed_step = propStep;

  if ( proposed_step >= retStep ) { // will cross into next volume next step
    moab::ErrorCode rval = DAG->next_vol(next_surf,vol,newvol);
    if(moab::MB_SUCCESS != rval) fludag_abort("g_fire","DAGMC failed in next_vol",rval);

    newRegion = DAG->index_by_handle(newvol);
    //      retStep = retStep; // path limited by geometry
    state.next_surface = next_surf; // no operation - but for clarity
    state.on_boundary = true;
    // history is preserved
  } else { // step less than the distance to surface or false step
    newRegion = oldRegion; // dont leave the current region
    retStep = proposed_step;    // physics limits step
    state.next_surface = state.prev_surface; // still hit the previous surface
    state.history.reset();       // reset the history
    state.on_boundary = false;     // cannot be on boundary
  }

  state.PrevRegion = newRegion; // particle will be moving to PrevRegion upon next entry.

  /*
  if ( mulbou_.lsense ) {
    // if a sensing step, then we will be fired again
    state.history.rollback_last_intersection();
    state.next_surf = state.prev_surf;
    state.on_boundary = false;
    retStep = mulbou_.tsense;
    // fluka takes care of reseting newreg outisde of ourtine
    //    proposed_step = mulbou_.tsense;
    if(debug)
      std::cout << "This is a sensing step " << mulbou_.tsense << std::endl;
  }
  */


  if(debug)  {
    std::cout << "Region on other side of surface is  = " << newRegion
              << ", Distance to next surf is " << retStep << std::endl;
  }

  // if lsense true, a false sensing step, facet not hit, surface not crossed
  // if sensing step takes is beyond volume doesnt matter, if takes in same volume
  // doesnt matter
  //  if(mulbou_.lsense) {
  //    state.history.rollback_last_intersection();
  //    if(debug) {
  //      std::cout << "This is a sensing step " << mulbou_.tsense << std::endl;
  //    }
  //    newRegion = -3;
  //  }

  // save all the state we need
  state.old_direction[0]=dir[0];
  state.old_direction[1]=dir[1];
  state.old_direction[2]=dir[2];

  // position
  state.old_position[0]=point[0];
  state.old_position[1]=point[1];
  state.old_position[2]=point[2];


  return;
}

// function to determine the angle between
double dot_product(moab::EntityHandle surface, double point[3], double direction[3])
{
  moab::EntityHandle rval;
  double normal[3] = {0,0,0};
  if(debug) {
    std::cout << "dot_product on surface " << surface << std::endl;
  }
  rval = DAG->get_angle(surface,point,normal);
  if(moab::MB_SUCCESS != rval) fludag_abort("dot_product","DAGMC failed in get_angle",rval);

  //may want to write our own inline vec angle function
  moab::CartVect v1(normal);
  moab::CartVect v2(direction);

  double angle = moab::angle(v1,v2);

  if(debug) {
    std::cout << "------ dot_product -------" << std::endl;
    std::cout << "angle between surface, " << surface << "with normal, " << v1 << std::endl;
    std::cout << " amd direction " << v2 << " is " << angle << std::endl;
  }

  return angle;
}

/* print state */
void print_state(particle_state &state)
{
  std::cout << "-------print_state-------" << std::endl;
  std::cout << "pos= " << state.old_position[0] << " " << state.old_position[1];
  std::cout << " "     << state.old_position[2] << std::endl;
  std::cout << "dir= " << state.old_direction[0] << " " << state.old_direction[1];
  std::cout << " "     << state.old_direction[2] << std::endl;

  std::cout << "next_surface: " << state.next_surface << " prev_surface: " << state.prev_surface << std::endl;
  std::cout << "prevRegion: " << state.PrevRegion  << std::endl;
  std::cout << "on_boundary: " << state.on_boundary << std::endl;
  return;
}

/* resets state */
void reset_state(particle_state &state)
{
  if(debug) {
    std::cout << "-----reset_state----" << std::endl;
    std::cout << "Resetting state....." << std::endl;
  }

  // set position to some large value
  state.old_position[0] = 9.e99;
  state.old_position[1] = 9.e99;
  state.old_position[2] = 9.e99;

  // direction set out of range of acceptable direction vectors
  state.old_direction[0] = -2.0;
  state.old_direction[1] = -2.0;
  state.old_direction[2] = -2.0;

  state.next_surface = 0;
  state.prev_surface = 0;
  state.PrevRegion = 0;

  state.on_boundary = false;

  state.history.reset();

  return;
}


/* testable wrapper for f_normal */
int normal (double& posx, double& posy, double& posz, double *norml, int& curRegion)
{
  int flagErr;
  int dummyReg;
  double dummyDirx, dummyDiry, dummyDirz;
  f_normal(posx, posy, posz, dummyDirx, dummyDiry, dummyDirz, norml, curRegion, dummyReg, flagErr);
  return flagErr;
}

/* given the particle position, direction, region return the normal to surface */
void f_normal(double& pSx, double& pSy, double& pSz,
              double& pVx, double& pVy, double& pVz,
              double* norml, const int& oldRegion,
              const int& newReg, int& flagErr)
{
  if(debug)
    std::cout << "============ NRMLWR =============" << std::endl;

  moab::EntityHandle OldReg = DAG -> entity_by_index(3,oldRegion); // entity handle
  double xyz[3] = {pSx,pSy,pSz}; //position vector
  double uvw[3] = {pVx,pVy,pVz}; //particl directoin
  int result; // particle is entering or leaving

  moab::ErrorCode rval = DAG->test_volume_boundary( OldReg, state.next_surface, xyz, uvw
                         ,result, &state.history);  // see if we are on boundary
  if(moab::MB_SUCCESS != rval) fludag_abort("f_normal","DAGMC failed in test_volume_boundary",rval);

  rval = DAG->get_angle(state.next_surface,xyz,norml);
  if(moab::MB_SUCCESS != rval) fludag_abort("f_normal","DAGMC failed in get_angle",rval);

  // result = 1 entering, 0 leaving
  if ( result == 0 ) { // vector should point towards OldReg
    norml[0] = norml[0]*-1.0;
    norml[1] = norml[1]*-1.0;
    norml[2] = norml[2]*-1.0;
  }

  if(debug)
    std::cout << "Normal: " << norml[0] << ", " << norml[1] << ", " << norml[2]  << std::endl;

  return;
}

/* does the position pos belong to region oldRegion */
inline bool check_vol( double pos[3], double dir[3], int oldRegion)
{
  int is_inside; // in volume or not
  // convert region id into entityhandle
  moab::EntityHandle volume = DAG->entity_by_index(3, oldRegion); // get the volume by index
  moab::ErrorCode rval = DAG->point_in_volume(volume, pos, is_inside,dir);

  // check for non error
  if(moab::MB_SUCCESS != rval) fludag_abort("check_vol","DAGMC failed in point_in_volume",rval);

  if ( is_inside == 1 ) // we are inside the cell tested
    return true;
  else
    return false;
}


/* testable wrapper function for f_look */
int look( double& posx, double& posy, double& posz, double* dir, int& oldRegion)
{
  int flagErr;
  int lattice_dummy;  // not used
  int nextRegion;
  f_look(posx, posy, posz, dir, oldRegion, lattice_dummy, nextRegion, flagErr, lattice_dummy);
  return nextRegion;
}


/* determine where a particle is, given position, direction etc */
void f_look(double& pSx, double& pSy, double& pSz,
            double* pV, const int& oldReg, const int& oldLttc,
            int& nextRegion, int& flagErr, int& newLttc)
{
  if(debug) {
    std::cout << "-------f_look------" << std::endl;
    std::cout << " stack size is " << flkstk_.npflka << std::endl;
    std::cout << "position is " << pSx << " " << pSy << " " << pSz << std::endl;
  }

  reset_state(state);
  // get the stack count
  //  state.stack_count = flkstk_.npflka;

  const double xyz[] = {pSx, pSy, pSz};       // location of the particle (xyz)
  const double dir[] = {pV[0],pV[1],pV[2]};

  int is_inside = 0;
  int num_vols = DAG->num_entities(3);  // number of volumes

  for (int i = 1 ; i <= num_vols ; i++) { // loop over all volumes
    moab::EntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
    // No ray history  - doesnt matter, only called for new source particles
    moab::ErrorCode rval = DAG->point_in_volume(volume, xyz, is_inside, dir);
    // check for non error
    if(moab::MB_SUCCESS != rval) fludag_abort("f_look","DAGMC failed in point_in_volume",rval);

    // we think that we are inside ,
    if ( is_inside == 1 ) { // we are inside the cell tested
      // test point in vol again, but reverse the direction
      double new_dir[3];
      new_dir[0] = -1.*dir[0];
      new_dir[1] = -1.*dir[1];
      new_dir[2] = -1.*dir[2];

      int second_test = 0 ;

      moab::ErrorCode rval = DAG->point_in_volume(volume, xyz, second_test, new_dir);
      // check for non error
      if(moab::MB_SUCCESS != rval) fludag_abort("f_look","DAGMC failed in point_in_volume",rval);

      if(is_inside != second_test) {
        if(debug) {
          std::cout << "First test inconclusive, doing the slow test" << std::endl;
        }

        rval =  DAG->point_in_volume_slow(volume, xyz, second_test);
        if(moab::MB_SUCCESS != rval) fludag_abort("f_look","DAGMC failed in point_in_volume_slow",rval);
        // firing rays in two opposing directions didnt work, the slow test didnt work
        if(second_test == 0 ) {
          continue;
        }
      }

      //WHEN WE ARE INSIDE A VOLUME, BOTH, nextRegion has to equal flagErr
      nextRegion = i;
      flagErr = nextRegion;

      if(debug) {
        std::cout << "region is " << nextRegion << " aka " << volume << std::endl;
      }
      return;
    } else if ( is_inside == -1 ) {
      std::cout << "We cannot be here" << std::endl;
      exit(0);
    }
  }  // end loop over all volumes

  // if are here then no volume has been found
  nextRegion = -33;
  flagErr = nextRegion;

  if(debug)
    std::cout <<  oldReg << " " <<  oldLttc << " " << nextRegion <<  " " << flagErr << " " << newLttc << std::endl;

  return;
}

/* If particle is lost calls this routine */
void f_lostlook(double& pSx, double& pSy, double& pSz,
                double* pV,  int& oldReg, int& oldLttc,
                int& nextRegion, int& flagErr, int& newLttc)
{
  if(debug) {
    std::cout << "-------f_lostlook------" << std::endl;
    std::cout << " stack size is " << flkstk_.npflka << std::endl;
    std::cout << "position is " << pSx << " " << pSy << " " << pSz << std::endl;
  }

  reset_state(state);
  // get the stack count
  //  state.stack_count = flkstk_.npflka;

  const double xyz[] = {pSx, pSy, pSz};       // location of the particle (xyz)
  const double dir[] = {pV[0],pV[1],pV[2]};

  int is_inside = 0;
  int num_vols = DAG->num_entities(3);  // number of volumes

  for (int i = 1 ; i <= num_vols ; i++) { // loop over all volumes
    moab::EntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
    // No ray history  - doesnt matter, only called for new source particles
    moab::ErrorCode rval = DAG->point_in_volume(volume, xyz, is_inside, dir);
    // check for non error
    if(moab::MB_SUCCESS != rval) fludag_abort("f_lostlook","DAGMC failed in point_in_volume",rval);

    // we think that we are inside ,
    if ( is_inside == 1 ) { // we are inside the cell tested
      // test point in vol again, but reverse the direction
      double new_dir[3];
      new_dir[0] = -1.*dir[0];
      new_dir[1] = -1.*dir[1];
      new_dir[2] = -1.*dir[2];

      int second_test = 0 ;

      moab::ErrorCode rval = DAG->point_in_volume(volume, xyz, second_test, new_dir);
      // check for non error
      if(moab::MB_SUCCESS != rval) fludag_abort("f_lostlook","DAGMC failed in point_in_volume",rval);

      if(is_inside != second_test) {
        if(debug) {
          std::cout << "First test inconclusive, doing the slow test" << std::endl;
        }

        rval =  DAG->point_in_volume_slow(volume, xyz, second_test);
        if(moab::MB_SUCCESS != rval) fludag_abort("f_lostlook","DAGMC failed in point_in_volume_slow",rval);
        // firing rays in two opposing directions didnt work, the slow test didnt work
        if(second_test == 0 ) {
          continue;
        }
      }

      //WHEN WE ARE INSIDE A VOLUME, BOTH, nextRegion has to equal flagErr
      nextRegion = i;
      flagErr = nextRegion;

      if(debug) {
        std::cout << "region is " << nextRegion << " aka " << volume << std::endl;
      }
      return;
    } else if ( is_inside == -1 ) {
      std::cout << "We cannot be here" << std::endl;
      exit(0);
    }
  }  // end loop over all volumes

  // if are here then no volume has been found
  nextRegion = num_vols + 1; // return nextRegion
  flagErr = nextRegion;

  if(debug)
    std::cout <<  oldReg << " " <<  oldLttc << " " << nextRegion <<  " " << flagErr << " " << newLttc << std::endl;

  return;
}

/* entering or leaving, if particle on boundary */
int boundary_test(moab::EntityHandle vol, double xyz[3], double uvw[3])
{
  int result;
  moab::ErrorCode rval = DAG->test_volume_boundary(vol,state.next_surface,
                         xyz,uvw, result,&state.history);  // see if we are on boundary
  // check for non error
  if(moab::MB_SUCCESS != rval) fludag_abort("boundary_test","DAGMC failed in test_volume_boundary",rval);

  return result;
}

/* Particle localisation when magnetic field tracking is on */
void lkmgwr(double& pSx, double& pSy, double& pSz,
            double* pV, const int& oldReg, const int& oldLttc,
            int& flagErr, int& newReg, int& newLttc)
{
  const double xyz[] = {pSx, pSy, pSz}; // location of the particle (xyz)
  int is_inside = 0; // logical inside or outside of volume
  int num_vols = DAG->num_entities(3); // number of volumes

  for (int i = 1 ; i <= num_vols ; i++) { // loop over all volumes
    moab::EntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
    // No ray history or ray direction.
    moab::ErrorCode rval = DAG->point_in_volume(volume, xyz, is_inside);

    // check for non error
    if(moab::MB_SUCCESS != rval) fludag_abort("lkmgwr","DAGMC failed in point_in_volume",rval);

    if ( is_inside == 1 ) { // we are inside the cell tested
      newReg = i;
      flagErr = i+1;
      if(debug) {
        std::cout << "point is in region = " << newReg << std::endl;
      }
      return;
    }
  }  // end loop over all volumes

  if(debug) {
    std::cout << "particle is nowhere!" << std::endl;
    std::cout << "point is not in any volume" << std::endl;
  }
  newReg = -33;
  return;
}

/* */
void f_lookdb(double& pSx, double& pSy, double& pSz,
              double* pV, const int& oldReg, const int& oldLttc,
              int& newReg, int& flagErr, int& newLttc)
{
  if(debug) {
    std::cout<<"============= F_LooKDB =============="<< std::endl;
  }

  //return region number and dummy variables
  newReg=0;
  newLttc=0;
  flagErr=-1;

  return;
}


/* terminates a history */
void f_g1rt(void)
{
  // reset all state
  reset_state(state);

  if(debug) {
    std::cout<<"============ F_G1RT ============="<<std::endl;
  }
  return;
}

/* Set DNEAR option if needed */
int f_idnr(const int & nreg, const int & mlat)
{
  if(debug)
    std::cout << "============= F_IDNR ==============" << std::endl;

  // returns 0 if user doesn't want Fluka to use DNEAR to compute the
  // step (the same effect is obtained with the GLOBAL (WHAT(3)=-1)
  // card in fluka input), returns 1 if user wants Fluka always to use DNEAR.
  return 0;
}

/* Wrapper for getting region name corresponding to given region number */
void rg2nwr(const int& mreg, const char* Vname)
{
  std::string vvname;

  if(debug) {
    std::cout << "============= RG2NWR ==============" << std::endl;
    std::cout << "mreg=" << mreg << std::endl;
  }

  region2name(mreg, vvname);
  Vname = vvname.c_str();

  if(debug) {
    std::cout << "reg2nmwr: Vname " << Vname<< std::endl;
  }

  return;
}

/* does nothing */
void rgrpwr(const int& flukaReg, const int& ptrLttc, int& g4Reg,
            int* indMother, int* repMother, int& depthFluka)
{
  if(debug)
    std::cout << "============= RGRPWR ==============" << std::endl;
  return;
}

/* returns magnetic field values in a given position */
void fldwr(const double& pX, const double& pY, const double& pZ,
           double& cosBx, double& cosBy, double& cosBz,
           double& Bmag, int& reg, int& idiscflag)

{
  if(debug)
    std::cout<<"================== MAGFLD ================="<<std::endl;
  return;
}

/* does nothing */
void flgfwr ( int& flkflg )
{
  if(debug)
    std::cout << "=======FLGFWR =======" << std::endl;
  return;
}

/* when lost try to fine the particle */
void lkfxwr(double& pSx, double& pSy, double& pSz,
            double* pV,  int& oldReg, int &oldLttc,
            int& newReg, int& flagErr, int& newLttc)
{
  if(debug)
    std::cout << "======= LKFXWR =======" << std::endl;

  f_lostlook(pSx,pSy,pSz,pV,oldReg,oldLttc,newReg,flagErr,newLttc);
  return;
}

/**************************************************************************************************/
/******                                End of FLUKA stubs                                  ********/
/**************************************************************************************************/

void fludag_abort(const char* function_name, const char* error_message, int error_code)
{
  std::stringstream ss;
  ss << std::string(error_message) << " " << error_code;
  std::string message = ss.str();
  const char* error_message_and_error = message.c_str();
  flabrt(function_name,error_message_and_error,strlen(function_name),strlen(error_message_and_error));
}

/**************************************************************************************************/
/******                                End of utility functions                            ********/
/**************************************************************************************************/

/* make the set of nuclides that are to be retained in full */
std::set<int> make_exception_set()
{
  std::set<int> nuc_exceptions;

  // Preserve FLUKA Entropy: Ref Fluka Manual pp. 318-323
  // Question about
  // Xenon (many named isotopes, useful only for detectors?)

  nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("H-1")));   // HYDROG-1
  nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("H-2")));   // DEUTERIU
  nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("H-3")));   // TRITIUM

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

  // All isotopes should be on the exception list, including the base isotope
  nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("H")));
  nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("He")));
  nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Li")));
  nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("B")));
  nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("Sr")));
  nuc_exceptions.insert(pyne::nucname::id(const_cast<char *>("I")));

  // Print out results
  if (debug) {
    std::cout << "Nucids of FLUKA exceptions" << std::endl;
    int i=1;
    for (std::set<int>::iterator ptr = nuc_exceptions.begin(); ptr != nuc_exceptions.end(); ++ptr) {
      std::cout << std::setw(10) << std::right << *ptr;
      if (i%5 == 0) {
        std::cout << std::endl;
      } else {
        std::cout << ", ";
      }
      i++;
    }
    std::cout << std::endl;
  }
  return nuc_exceptions;
}

// FluDAG Material Card  Functions
void fludag_write_ididx(std::string ididx)
{
  int num_vols = DAG->num_entities(3);
  std::ofstream file;
  file.open(ididx.c_str());
  for ( int i = 1 ; i <= num_vols ; i++ ) {
    int cellid = DAG->id_by_index( 3, i );
    moab::EntityHandle entity = DAG->entity_by_index( 3, i );
    file << i << " " << cellid << " " << entity << std::endl;
  }
  file.close();
}

// FluDAG Material Card  Functions
void fludag_write(std::string matfile, std::string lfname)
{


  // Use DAG to read and count the volumes.
  std::map<int, std::string> map_name;
  if (0 == DAG->num_entities(3) ) {
    std::cout << "Error: there are no volumes in this geometry!" << std::endl;
    return;
  }

  // get the pyne materials and tallies
  UWUW workflow_data = UWUW(matfile);

  std::list<pyne::Material> pyne_list;
  std::map<std::string, pyne::Material> pyne_map;
  pyne_map = workflow_data.material_library;

  // ASSIGNMA Cards
  std::ostringstream astr;
  fludagwrite_assignma(astr, pyne_map, map_name);

  // MATERIAL Cards
  pyne::NUC_DATA_PATH = workflow_data.full_filepath; // for atomic data

  // write COMPOUND CARDS
  std::ostringstream mstr;
  fludag_all_materials(mstr, pyne_map);

  // Write all the streams to the input file
  std::string header = "*...+....1....+....2....+....3....+....4....+....5....+....6....+....7...";
  std::ofstream lcadfile (lfname.c_str());
  lcadfile << header << std::endl;
  lcadfile << astr.str();
  lcadfile << header << std::endl;
  lcadfile << mstr.str();

  lcadfile << header << std::endl;
  lcadfile << "* UW**2 tallies" << std::endl;
  mstr.str("");
  fludag_all_tallies(mstr,workflow_data.tally_library);
  lcadfile << mstr.str();

  // all done
  lcadfile.close();
}

//---------------------------------------------------------------------------//
// fludagwrite_assignma
//---------------------------------------------------------------------------//
// Put the ASSIGNMAt statements in the output ostringstream
void fludagwrite_assignma(std::ostringstream& ostr,
                          std::map<std::string, pyne::Material> pyne_map,
                          std::map<int, std::string> map_name)
{

  // get the material and density props
  std::map<moab::EntityHandle,std::vector<std::string> > material_assignments = get_property_assignments("mat",3,":/");
  std::map<moab::EntityHandle,std::vector<std::string> > density_assignments = get_property_assignments("rho",3,":/");

  pyne::Material material;

  std::vector<std::string> material_props,density_props;

  // loop over all volumes
  for (unsigned int vol_i = 1 ; vol_i <= DAG->num_entities(3) ; vol_i++) {
    int cellid = DAG->id_by_index( 3, vol_i );
    moab::EntityHandle entity = DAG->entity_by_index( 3, vol_i );

    material_props = material_assignments[entity];
    density_props = density_assignments[entity];

    if( material_props.size() > 1 ) {
      std::cout << "more than one material for volume with id " << cellid << std::endl;
      std::cout << cellid << " has the following material assignments" << std::endl;
      for ( int j = 0 ; j < material_props.size() ; j++ ) {
        std::cout << material_props[j] << std::endl;
      }
      std::cout << "Please check your material assignments " << cellid << std::endl;
      exit(EXIT_FAILURE);
    }

    if(density_props.size() > 1) {
      std::cout << "More than one density specified for " << cellid <<std::endl;
      std::cout << cellid << " has the following density assignments" << std::endl;
      for ( int j = 0 ; j < density_props.size() ; j++ ) {
        std::cout << density_props[j] << std::endl;
      }
      std::cout << "Please check your density assignments " << cellid << std::endl;
      exit(EXIT_FAILURE);
    }

    std::string grp_name = "";
    if (!density_props[0].empty()) {
      grp_name = "mat:"+material_props[0]+"/rho:"+density_props[0];
    } else {
      grp_name = "mat:"+material_props[0];
    }

    std::string fluka_name = "";

    // not graveyard or vacuum or implicit compliment
    if (grp_name.find("Graveyard") == std::string::npos && grp_name.find("Vacuum") == std::string::npos
        && !(DAG->is_implicit_complement(entity)) ) {
      material = pyne_map[grp_name];
      fluka_name = material.metadata["fluka_name"].asString();
    } else if (grp_name.find("Graveyard") != std::string::npos ||
               grp_name.find("graveyard") != std::string::npos ) {

      fluka_name = "BLCKHOLE";
    } else if (grp_name.find("Vacuum") != std::string::npos) {
      fluka_name = "VACUUM";
    }  else if (  DAG->is_implicit_complement(entity) ) {
      fluka_name = "VACUUM";
    }

    // The fluka name has been found, create the card
    ostr << std::setw(10) << std::left  << "ASSIGNMA ";
    ostr << std::setw(10) << std::right << fluka_name;
    ostr << std::setprecision(0) << std::fixed << std::showpoint
         << std::setw(10) << std::right << (float)vol_i << std::endl;

  }   // End loop through vol_i
  std::cout << std::endl;

}  // end fludagwrite_assignma


// Get tally cards for all tallies in the problem
void fludag_all_tallies(std::ostringstream& mstr, std::map<std::string,pyne::Tally> tally_map)
{
  int start_unit = 21; // starting unit number for tallies

  std::map<std::string,pyne::Tally>::iterator it;

  // generate number of tally/particle pairs
  std::list<std::string> tally_parts;
  std::string tally_id;
  for ( it = tally_map.begin() ; it != tally_map.end() ; ++it ) {
    tally_id = (it->second).tally_type+"/"+(it->second).particle_name;
    if( std::count(tally_parts.begin(),tally_parts.end(),tally_id) == 0 ) {
      tally_parts.insert(tally_parts.end(),tally_id);
    }
  }

  // loop over tallies in map
  for ( it = tally_map.begin() ; it != tally_map.end() ; ++it ) {
    pyne::Tally tally = (it->second);
    // pyne tallies are by id, FluDAG is by index, need to convert
    moab::EntityHandle vol_eh = DAG->entity_by_id(3,tally.entity_id);
    // volume index
    int vol_idx = DAG->index_by_handle(vol_eh);
    // recast tally to index, use entity_name for setting volume

    moab::ErrorCode rval = DAG->measure_volume(vol_eh,tally.entity_size);

    std::stringstream ss;
    ss << vol_idx;
    ss << ".";
    tally.entity_name = ss.str();

    std::string tally_id = tally.tally_type+"/"+tally.particle_name;

    std::list<std::string>::iterator iter = std::find (tally_parts.begin(), tally_parts.end(), tally_id);

    int unit_number = std::distance(tally_parts.begin(), iter) + start_unit;

    ss.str(std::string());
    ss << "-";
    ss << unit_number;


    mstr << tally.fluka(ss.str()) << std::endl;
  }

  return;
}

//---------------------------------------------------------------------------//
// fludag_all_materials
//---------------------------------------------------------------------------//
// Get material cards for all materials in the problem, both elemental and compounds
void fludag_all_materials(std::ostringstream& mstr, std::map<std::string,pyne::Material> pyne_map)
{
  std::set<int> exception_set = make_exception_set();

  std::map<int, std::string> map_nucid_fname;

  pyne::Material unique = pyne::Material();

  // loop over all materials, summing
  std::map<std::string, pyne::Material>::iterator nuc;
  for ( nuc = pyne_map.begin(); nuc != pyne_map.end(); ++nuc) {
    unique = unique + (nuc->second);
  }
  // now collapse elements
  unique = unique.collapse_elements(exception_set);

  // remove those that are no longer needed due to
  // compound card inclusions
  // now write out material card & compound card for each compound
  for ( nuc = pyne_map.begin() ; nuc != pyne_map.end(); ++nuc) {
    pyne::Material compound = (nuc->second).collapse_elements(exception_set);
    // if only one element in comp, then we can remove the one that exists
    // in the unique material
    if ( compound.comp.size() == 1 ) {
      pyne::comp_iter nuc = compound.comp.begin();
      std::set<int> nuc2del;
      nuc2del.insert(nuc->first);
      // remove the nuclide from the unique list
      unique = unique.del_mat(nuc2del);
    }
  }

  //del_mat

  // number of required material cards due to calls
  int num_mat = unique.comp.size();

  // write out material card for each one
  int i = ID_START;
  pyne::comp_map::iterator element;
  std::string mat_line;
  for ( element = unique.comp.begin() ; element != unique.comp.end() ; ++element) {
    int nuc_id = element->first; // get the nuc id
    pyne::comp_map nucvec;
    nucvec[nuc_id] = 100.0; // create temp nucvec
    pyne::Material element_tmp = pyne::Material(nucvec); // create temp material
    mat_line = element_tmp.fluka(i);
    if (mat_line.length() != 0) {
      i++;
      mstr << mat_line;
    }
  }

  // now write out material card & compound card for each compound
  std::string compound_string;
  for ( nuc = pyne_map.begin() ; nuc != pyne_map.end(); ++nuc) {
    pyne::Material compound = (nuc->second).collapse_elements(exception_set);
    compound_string = compound.fluka(i);
    if ( compound_string.length() != 0 ) {
      i++;
      mstr << compound_string;
    }
  }
}

// region2name
void region2name(int volindex, std::string &vname )  // file with cell/surface cards
{
  std::stringstream ss;
  ss << volindex;
  ss << ".";
  vname = ss.str();
}

// get all property in all volumes
std::map<moab::EntityHandle,std::vector<std::string> > get_property_assignments(std::string property,
    int dimension, std::string delimiters)
{

  std::map<moab::EntityHandle,std::vector<std::string> > prop_map;

  std::vector< std::string > mcnp5_keywords;
  std::map< std::string, std::string > mcnp5_keyword_synonyms;

  // populate keywords
  mcnp5_keywords.push_back( "mat" );
  mcnp5_keywords.push_back( "rho" );
  mcnp5_keywords.push_back( "tally" );

  // get initial sizes
  int num_entities = DAG->num_entities( dimension );

  // parse data from geometry
  moab::ErrorCode rval = DAG->parse_properties( mcnp5_keywords, mcnp5_keyword_synonyms,delimiters.c_str());

  if (moab::MB_SUCCESS != rval) {
    std::cout << "DAGMC failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }


  // loop over all cells
  for( int i = 1; i <= num_entities; ++i ) {
    // get cellid
    moab::EntityHandle entity = DAG->entity_by_index( dimension, i );

    std::vector<std::string> properties;
    std::vector<std::string> tmp_properties;


    // get the group contents
    if( DAG->has_prop( entity, property ) ) {
      rval = DAG->prop_values(entity,property,tmp_properties);
      properties.push_back(tmp_properties[0]);
    } else
      properties.push_back("");

    prop_map[entity]=properties;

  }

  return prop_map;
}
